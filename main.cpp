#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "DSUeval.h"


extern "C" {
  #include "DSUeval_cmdline.h"
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
  if (args_info.print_full_flag || args_info.energy_file_given) args_info.print_energy_flag = 1;

  // code
    // DSUeval
  DSU dsu(stdin, args_info.noLP_flag, args_info.shift_flag, args_info.time_max_arg, args_info.number_lm_arg);
  Opt opt(args_info.noLP_flag, args_info.shift_flag, !args_info.noSaddle_flag, args_info.floodMax_arg, args_info.floodHeight_arg, args_info.debug_flag, args_info.depth_arg, args_info.num_threshold_arg, args_info.outer_flag, args_info.cluster_repre_arg, !args_info.cluster_fsaddle_flag);

  // adjust args_info:
  if (args_info.cluster_Kmax_arg>=dsu.Size()) {
    args_info.cluster_off_flag = 1;
  }

  dsu.Cluster(opt, args_info.cluster_Kmax_arg, args_info.cluster_off_flag);
  fprintf(stderr, "computation of saddles took %.2f secs.\n", (clock()-time)/(double)CLOCKS_PER_SEC);

  if (args_info.just_ub_flag) {
    printf("%.2f %d\n", (clock()-time)/(double)CLOCKS_PER_SEC, dsu.Size()); time = clock();

  } else {
    time = clock();
      // LinkCP
    dsu.LinkCPLM(opt, args_info.debug_flag);
    fprintf(stderr, "computing LM-edges took %.2f secs.\n", (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();
    if (!args_info.noSaddle_flag) {
      dsu.LinkCPsaddle(opt, args_info.debug_flag);
      fprintf(stderr, "computing saddle-edges took %.2f secs.\n", (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();
    }
    dsu.SortFix();
    dsu.PrintDot(args_info.name_dot_arg, args_info.dot_flag, args_info.print_graph_flag, args_info.name_graph_arg, args_info.tree_visualise_flag);

    fprintf(stderr, "printing dot took %.2f secs.\n", (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();

      // connect comps
    if (args_info.components_flag) {
      dsu.PrintComps();
      dsu.ConnectComps(args_info.depth_arg, args_info.debug_flag);
      //dsu.PrintDot(args_info.name_dot_arg, args_info.dot_flag, args_info.print_graph_flag, args_info.name_graph_arg, args_info.tree_visualise_flag);
      dsu.PrintDot(args_info.name_dot_arg, args_info.dot_flag, args_info.print_graph_flag, args_info.name_graph_arg, args_info.tree_visualise_flag);
      fprintf(stderr, "connnecting components took %.2f secs.\n", (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();
    }
    /*// just debug
    dsu.PrintComps();
    dsu.PrintLinkCP(false);
  */


    // real output:
    dsu.PrintLM(stdout);

    // barriers-like output
    if (args_info.barr_file_given) {
      FILE *file_h;
      file_h = fopen(args_info.barr_file_arg, "w");
      dsu.PrintBarr(file_h);
      fclose(file_h);
    }

    //dsu.PrintSaddles(stderr);
    //dsu.PrintComps(stderr, true);

    // print energy matrix
    if (args_info.print_energy_flag) {
      dsu.PrintMatrix(args_info.energy_file_arg, args_info.print_full_flag, 'E');
    }

    // print dist matrix
    if (args_info.dist_file_given) {
      dsu.PrintMatrix(args_info.dist_file_arg, args_info.print_full_flag, 'D');
    }

    // print rates matrix
    if (args_info.rates_file_given) {
      dsu.PrintRates(args_info.rates_file_arg, args_info.print_full_flag, args_info.rates_temp_arg, args_info.rates_mode_arg[0]);
    }

    // print graph distance matrix
    if (args_info.gdist_file_given) {
      dsu.PrintMatrix(args_info.gdist_file_arg, args_info.print_full_flag, 'G');
    }

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
      FILE *file_h;
      file_h = fopen(args_info.energy_heights_arg, "w");
      bool full = true;
      dsu.EHeights(file_h, full);
      fclose(file_h);
    }

    // evaluation
    if (args_info.energy_rank_given) {
      FILE *file_h;
      file_h = fopen(args_info.energy_rank_arg, "w");
      dsu.ERank(file_h, false, true);
      fclose(file_h);
    }

    // evaluation
    if (args_info.energy_barrier_given) {
      FILE *file_h;
      file_h = fopen(args_info.energy_barrier_arg, "w");
      dsu.ERank(file_h, true);
      fclose(file_h);
    }
  }

  fprintf(stderr, "DSUeval exitting succesfully!\n");

  cmdline_parser_free(&args_info);
  fprintf(stderr, "rest took %.2f secs.\n", (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();
}
