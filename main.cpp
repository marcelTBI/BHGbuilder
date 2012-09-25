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
  if (args_info.num_threshold_arg <= 0) args_info.hd_threshold_arg = INT_MAX;

  // code
    // DSUeval
  DSU dsu(stdin);
  dsu.CreateList(args_info.hd_threshold_arg, args_info.debug_flag);
  dsu.ComputeUB(args_info.depth_arg, args_info.num_threshold_arg, args_info.outer_flag, args_info.debug_flag);
  dsu.PrintUBoutput();

    // LinkCP
  dsu.LinkCP(args_info.shift_flag, args_info.noLP_flag, args_info.debug_flag);
  dsu.PrintDot(args_info.name_dot_arg, args_info.dot_flag, args_info.print_graph_flag, args_info.name_graph_arg, args_info.tree_visualise_flag);
  if (args_info.print_energy_flag) dsu.PrintMatrix(args_info.energy_file_arg);

  cmdline_parser_free(&args_info);
}
