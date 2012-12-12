#include <stdio.h>
#include <stdlib.h>

#include "DSUeval.h"

extern "C" {
  #include "DSUeval_cmdline.h"
}

using namespace std;

int main(int argc, char **argv)
{
  // parse arguments
  gengetopt_args_info args_info;
  if (cmdline_parser(argc, argv, &args_info) != 0) {
    fprintf(stderr, "Argument parsing problem.\n");
    exit(EXIT_FAILURE);
  }

  //adjust args_info
  if (args_info.hd_threshold_arg <= 0) args_info.hd_threshold_arg = INT_MAX;
  if (args_info.num_threshold_arg <= 0) args_info.num_threshold_arg = INT_MAX;

  // code
    // DSUeval
  DSU dsu(stdin, args_info.noLP_flag, args_info.shift_flag);
  dsu.CreateList(args_info.hd_threshold_arg, args_info.debug_flag);
  dsu.ComputeUB(args_info.depth_arg, args_info.num_threshold_arg, args_info.outer_flag, args_info.noLP_flag, args_info.shift_flag, args_info.debug_flag);
  dsu.PrintUBoutput(stderr);

    // LinkCP
  Opt opt(args_info.noLP_flag, args_info.shift_flag, !args_info.noSaddle_flag, args_info.floodMax_arg, args_info.floodHeight_arg);
  dsu.LinkCP(opt, args_info.debug_flag);
  dsu.PrintDot(args_info.name_dot_arg, args_info.dot_flag, args_info.print_graph_flag, args_info.name_graph_arg, args_info.tree_visualise_flag);

    // connect comps
  if (args_info.components_flag) {
    //dsu.PrintComps();
    dsu.ConnectComps(args_info.depth_arg, args_info.debug_flag);
    //dsu.PrintDot(args_info.name_dot_arg, args_info.dot_flag, args_info.print_graph_flag, args_info.name_graph_arg, args_info.tree_visualise_flag);
  }
  /*// just debug
  dsu.PrintComps();
  dsu.PrintLinkCP(false);
*/
  // real output:
  dsu.PrintLM(stdout);
  dsu.PrintSaddles(stderr);
  dsu.PrintComps(stderr, true);

  dsu.PrintDot(args_info.name_dot_arg, args_info.dot_flag, args_info.print_graph_flag, args_info.name_graph_arg, args_info.tree_visualise_flag);


  // print energy matrix
  if (args_info.print_energy_flag) dsu.PrintMatrix(args_info.energy_file_arg);

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

  // temporary
  if (args_info.energy_heights_flag) {
    FILE *file_h;
    file_h = fopen("energy_heights.txt", "w");
    bool full = true;
    dsu.EHeights(file_h, full);
    fclose(file_h);
  }

  fprintf(stderr, "DSUeval exitting succesfully!\n");

  cmdline_parser_free(&args_info);
}
