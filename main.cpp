#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "BHGbuilder.h"
#include "RateGraph.h"


extern "C" {
  #include "BHGbuilder_cmdline.h"
}

using namespace std;

int main(int argc, char **argv)
{
  clock_t time = clock();
  // parse arguments
  gengetopt_args_info args_info;
  if (cmdline_parser(argc, argv, &args_info) != 0) {
    fprintf(stderr, "Argument parsing problem.\n");
    exit(EXIT_FAILURE);
  }

  //adjust args_info
  if (args_info.hd_threshold_arg <= 0) args_info.hd_threshold_arg = INT_MAX;
  if (args_info.num_threshold_arg <= 0) args_info.num_threshold_arg = INT_MAX;
  if (args_info.cluster_repre_arg <= 0.0) args_info.cluster_repre_arg = 0.0;
  if (args_info.cluster_repre_arg > 1.0) args_info.cluster_repre_arg = 1.0;

  // print-all
  if (args_info.print_all_flag) {
    args_info.barr_file_given = 1; free(args_info.barr_file_arg); args_info.barr_file_arg = (char*) malloc(50*sizeof(char));  strcpy(args_info.barr_file_arg, "barrier.lm");
    args_info.energy_file_given = 1; free(args_info.energy_file_arg); args_info.energy_file_arg = (char*) malloc(50*sizeof(char));  strcpy(args_info.energy_file_arg, "enbarr.eb");
    args_info.dist_file_given = 1; free(args_info.dist_file_arg); args_info.dist_file_arg = (char*) malloc(50*sizeof(char));  strcpy(args_info.dist_file_arg, "bpdistances.dist");
    args_info.gdist_file_given = 1; free(args_info.gdist_file_arg); args_info.gdist_file_arg = (char*) malloc(50*sizeof(char));  strcpy(args_info.gdist_file_arg, "graphdistances.gdist");
    args_info.rates_file_given = 1; free(args_info.rates_file_arg); args_info.rates_file_arg = (char*) malloc(50*sizeof(char));  strcpy(args_info.rates_file_arg, "rates.rat");
  }

  // code
    // BHGbuilder
  Opt opt(args_info);
  DSU dsu(stdin, args_info.noLP_flag, args_info.shift_flag, args_info.pseudoknots_flag, args_info.time_max_arg, args_info.number_lm_arg, args_info.just_read_flag, args_info.debug_flag);
  fprintf(stderr, "reading input took %.2f secs.\n", (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();

  if (!args_info.just_read_flag) {
    dsu.Cluster(opt, args_info.cluster_Kmax_arg);
    fprintf(stderr, "computation of saddles took %.2f secs.\n", (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();
  }

  if (args_info.just_ub_flag) {
    printf("%.2f %d\n", (clock()-time)/(double)CLOCKS_PER_SEC, dsu.Size()); time = clock();

  } else {
    if (!args_info.just_read_flag) {
      time = clock();
        // LinkCP
      dsu.LinkCPLM(opt, args_info.debug_flag);
      fprintf(stderr, "computing LM-edges took %.2f secs.\n", (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();
      if (!args_info.noSaddle_flag && args_info.floodMax_arg>0) {
        dsu.LinkCPsaddle(opt, args_info.debug_flag);
        fprintf(stderr, "computing saddle-edges took %.2f secs.\n", (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();
      }
      dsu.SortFix();

      fprintf(stderr, "sorting took %.2f secs.\n", (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();

        // connect comps
      if (args_info.components_flag) {
        //dsu.PrintComps(stdout);
        dsu.ConnectComps(args_info.depth_arg, args_info.debug_flag);
        //dsu.PrintDot(args_info.name_dot_arg, args_info.dot_flag, args_info.print_graph_flag, args_info.name_graph_arg, args_info.tree_visualise_flag);
        fprintf(stderr, "connnecting components took %.2f secs.\n", (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();
      }
    }

    // reduce!
    if (args_info.reduce_given) {
      dsu.Reduce(args_info.reduce_arg, args_info.filter_file_given?args_info.filter_file_arg:NULL);
      fprintf(stderr, "reduction took %.2f secs.\n", (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();
    }

    // print dot
    dsu.PrintDot(args_info.dot_file_arg, args_info.dot_flag, args_info.graph_file_given, args_info.graph_file_arg, args_info.tree_visualise_flag, args_info.dot_energies_flag);
    fprintf(stderr, "printing dot took %.2f secs.\n", (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();


    if (!args_info.quiet_flag) {
      // real output:
      dsu.PrintLM(stdout);
      dsu.PrintSaddles(stdout, true, args_info.no_new_flag);
    }

    // barriers-like output
    if (args_info.barr_file_given) {
      FILE *file_h;
      file_h = fopen(args_info.barr_file_arg, "w");
      dsu.PrintBarr(file_h);
      fclose(file_h);
      fprintf(stderr, "printing barrier output took %.2f secs.\n", (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();
    }

    //dsu.PrintSaddles(stderr);
    //dsu.PrintComps(stderr, true);

    dsu.SetkT(args_info.rates_temp_arg);

    if (args_info.bulk_path_given) {
      vector <int> indices_path;
      indices_path = dsu.GetNumbers(args_info.bulk_path_arg);
      for (unsigned int i=1; i<indices_path.size(); i++) {
        dsu.GetPath(indices_path[i-1], indices_path[i], args_info.depth_arg, args_info.optimal_path_arg[0], args_info.saddle_path_given?en_fltoi(args_info.saddle_path_arg):10000, args_info.time_path_arg, i);
      }
    }

    //print optimal path
    for (int i=0; i<(int)args_info.get_path_given; i++) {
      int a, b;
      if (sscanf(args_info.get_path_arg[i], "%d=%d", &a, &b)!=2) {
        fprintf(stderr, "WARNING: wrong use of --get-path option (%s)\n", args_info.get_path_arg[i]);
      } else {
        if (a<=0 || b<=0) {
          fprintf(stderr, "WARNING: non-positive number in --get-path (%s)\n", args_info.get_path_arg[i]);
        } else {
          a--;
          b--;
          if (a>=dsu.Size() || b>=dsu.Size()) {
            fprintf(stderr, "WARNING: visualisation number(s) exceeds number of minima (%d) (%s)\n", dsu.Size(), args_info.get_path_arg[i]);
          } else dsu.GetPath(a, b, args_info.depth_arg, args_info.optimal_path_arg[0], args_info.saddle_path_given?en_fltoi(args_info.saddle_path_arg):10000, args_info.time_path_arg);
        }
      }
    }

    //#######  this part is influenced by filter-file:

    // print energy matrix
    if (args_info.energy_file_given) {
      dsu.PrintMatrix(args_info.energy_file_arg, args_info.print_full_flag, args_info.filter_file_given?args_info.filter_file_arg:NULL, 'E');
    }

    // print dist matrix
    if (args_info.dist_file_given) {
      dsu.PrintMatrix(args_info.dist_file_arg, args_info.print_full_flag, args_info.filter_file_given?args_info.filter_file_arg:NULL, 'D');
    }

    // print bdist matrix
    if (args_info.bdist_file_given) {
      dsu.PrintMatrix(args_info.bdist_file_arg, args_info.print_full_flag, args_info.filter_file_given?args_info.filter_file_arg:NULL, 'B');
    }

    // print graph distance matrix
    if (args_info.gdist_file_given) {
      dsu.PrintMatrix(args_info.gdist_file_arg, args_info.print_full_flag, args_info.filter_file_given?args_info.filter_file_arg:NULL, 'G');
    }

    // print pknot-type matrix
    if (args_info.ptype_file_given) {
      dsu.PrintMatrix(args_info.ptype_file_arg, args_info.print_full_flag, args_info.filter_file_given?args_info.filter_file_arg:NULL, 'P');
    }

    fprintf(stderr, "printing matrices took %.2f secs.\n", (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();


    // Shur removal
    if (!args_info.filter_given && args_info.filter_file_given) {
      if (!args_info.filter_arg) free(args_info.filter_arg);
      args_info.filter_arg = (char*)malloc(sizeof(char)*(strlen(args_info.filter_file_arg)+1));
      strcpy(args_info.filter_arg, args_info.filter_file_arg);
      args_info.filter_given = 1;
    }
    if (args_info.rates_file_given || args_info.prob_param_given || args_info.fret_given) {
      RateGraph rg(dsu, args_info.rates_temp_arg, args_info.rates_fullpath_flag?args_info.depth_arg:0, args_info.minimal_rate_arg, args_info.ordering_arg[0]);

      // read filter
      if (args_info.filter_given) {
        args_info.max_arg = max(rg.ReadFilter(args_info.filter_arg), args_info.max_arg);
      }
      if (args_info.max_arg == -1 || args_info.max_arg > rg.Size()) args_info.max_arg = rg.Size();
      int x= 0;

      // start actual removing: remove one by one using the sparsity:
      if (rg.Size()-args_info.max_arg > 0) {
        int queue = rg.ConstructQueue(args_info.ordering_arg[0], rg.Size()-args_info.max_arg, args_info.leave_transitions_flag);
        fprintf(stderr, "creating graph took %.2f secs.\n", (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();
        x = rg.RemoveX(rg.Size()-args_info.max_arg, args_info.fraction_arg, !args_info.nreeval_flag, args_info.schur_maximal_arg, args_info.minimal_rate_arg);
        fprintf(stderr, "removal of %d lm took %.2f secs. (one by one removal)\n", x, (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();
      }

      // bulk remove:
      if (rg.Size()-args_info.max_arg-x > 0) {
        x = rg.RemoveShur(rg.Size()-args_info.max_arg-x, args_info.Schur_step_arg, args_info.minimal_rate_arg);
        fprintf(stderr, "removal of %d lm took %.2f secs. (bulk matrix removal)\n", x, (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();
      }

      if (args_info.rates_file_given) {
        // print rates
        FILE *file = fopen(args_info.rates_file_arg, "w");
        if (file) {
          fprintf(stderr, "printing rates ... ");
          rg.PrintRates(file);
          fprintf(stderr, "into \"%s\" took %.2f secs.\n", args_info.rates_file_arg, (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();
          fclose(file);
        }

        // now output the list of stuff that we haven't removed
        char outrates[500];
        strcpy(outrates, args_info.rates_file_arg);
        int len = strlen(outrates);
        outrates[len] = 'O';
        outrates[len+1] = '\0';
        rg.PrintOutput(outrates);
        //rg.PrintDot(args_info.dot_file_arg, args_info.graph_file_given);
      }

      // now output the prob matrix:
      if (args_info.prob_param_given) {
        char filename[500];
        sprintf(filename, "prob_matrix%.0f.tsv", args_info.prob_param_arg);
        fprintf(stderr, "exponentiating rate matrix ... ");
        rg.PrintProb(args_info.prob_param_arg, filename);
        fprintf(stderr, "into \"%s\" took %.2f secs.\n", filename, (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();
      }

      // now output all the FRET data:
      if (args_info.fret_given) {
        fprintf(stderr, "exponentiating rate matrix ... ");
        rg.CreateProb(args_info.prob_param_arg);

        // first file:
        char filename[500];
        sprintf(filename, "idstr_%s.tsv", args_info.fret_arg);
        dsu.PrintFRET1(filename);

        // second file:
        sprintf(filename, "iniprob_%s.tsv", args_info.fret_arg);
        dsu.PrintFRET2(filename);

        // third file:
        sprintf(filename, "tranprob_%s.tsv", args_info.fret_arg);
        rg.PrintProb(args_info.prob_param_arg, filename);

        fprintf(stderr, "creating FRET files took %.2f secs.\n", (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();
      }

    }

    //##########  end of influence

    /*// print rates matrix
    if (args_info.rates_file_given) {
      for (int i=0; i<max(1, (int)args_info.rates_mode_given); i++) {
        char filename[strlen(args_info.rates_file_arg)+2];
        strcpy(filename, args_info.rates_file_arg);
        filename[strlen(args_info.rates_file_arg)]=args_info.rates_mode_arg[i][0];
        filename[strlen(args_info.rates_file_arg)+1]='\0';
        //fprintf(stderr, filename);
        dsu.PrintRates(filename, args_info.print_full_flag, args_info.rates_temp_arg, args_info.rates_mode_arg[i][0]);
      }
      fprintf(stderr, "printing rates(%d) took %.2f secs.\n", max(1, (int)args_info.rates_mode_given), (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();
    }*/

    // visualisation
    for (unsigned int i=0; i<args_info.visualise_given; i++) {
      int a, b;
      if (sscanf(args_info.visualise_arg[i], "%d=%d", &a, &b)!=2) {
        fprintf(stderr, "WARNING: wrong use of visualisation option (%s)\n", args_info.visualise_arg[i]);
      } else {
        if (a<=0 || b<=0) {
          fprintf(stderr, "WARNING: non-positive number in visualisation (%s)\n", args_info.visualise_arg[i]);
        } else {
          a--;
          b--;
          if (a>=dsu.Size() || b>=dsu.Size()) {
            fprintf(stderr, "WARNING: visualisation number(s) exceeds number of minima (%d) (%s)\n", dsu.Size(), args_info.visualise_arg[i]);
          } else dsu.VisPath(a, b, !args_info.vis_dist_flag, args_info.vis_length_arg, args_info.dot_flag, args_info.debug_flag);
        }
      }
    }

    // evaluation
    if (args_info.energy_heights_given) {
      fprintf(stderr, "Printing energy heights...");
      FILE *file_h;
      file_h = fopen(args_info.energy_heights_arg, "w");
      bool full = true;
      bool only_norm = true;
      dsu.EHeights(file_h, full, only_norm);
      fclose(file_h);
      fprintf(stderr, "done, it took %.2f secs.\n", (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();
    }

    // evaluation
    if (args_info.energy_rank_given) {
      fprintf(stderr, "Printing energy ranks...");
      FILE *file_h;
      file_h = fopen(args_info.energy_rank_arg, "w");
      dsu.ERank(file_h, false, true);
      fclose(file_h);
      fprintf(stderr, "done, it took %.2f secs.\n", (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();
    }

    // evaluation
    if (args_info.energy_barrier_given) {
      fprintf(stderr, "Printing energy barriers...");
      FILE *file_h;
      file_h = fopen(args_info.energy_barrier_arg, "w");
      dsu.ERank(file_h, true, true);
      fclose(file_h);
      fprintf(stderr, "done, it took %.2f secs.\n", (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();
    }

    // graph analyze:
    if (args_info.analyze_graph_flag) {
      fprintf(stderr, "Analysing the graph...");
      FILE *histo;
      histo = fopen("histogram.txt", "w");
      dsu.Histo(histo);
      fclose(histo);
      fprintf(stderr, "done, it took %.2f secs.\n", (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();
    }
  }

  fprintf(stderr, "BHGbuilder exitting succesfully!\n");

  cmdline_parser_free(&args_info);
  fprintf(stderr, "rest took %.2f secs.\n", (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();
}
