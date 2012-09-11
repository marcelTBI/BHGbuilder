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
  dsu.PrintDot("graph.dot", true, true, "graph.eps", args_info.tree_visualise_flag);
  dsu.PrintMatrix("graph.en");

  cmdline_parser_free(&args_info);
}
