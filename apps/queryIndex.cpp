
#include <ConwayBromageDatastructure.h>
#include <getopt.h>

#include <fstream>
#include <iostream>

void usage(std::string name) {
  std::cerr << "Query a saved Conway-Bromage structure with a list of kmers.\n"
            << "Usage: " << name << " [OPTIONS] <index> <kmers>\n"
            << "\n"
            << "[ARGUMENTS]\n"
            << "    <index>  CBD containing file\n"
            << "    <kmers>  list of query k-mers (k-1 w.r.t to the ones used "
               "to build the graph)\n"
            << "\n"
            << "[OPTIONS]\n"
            << "    -s --successor  Do successor queries instead of contains\n"
            << "    -h --help       Show this message" << std::endl;
};

struct CliArgs {
  bool successor;
  std::string index_filename;
  std::string kmers_filename;
};

int parse_args(int argc, char* argv[], CliArgs& args) {
  // Arg option parsing
  int c;
  static struct option long_opts[] = {{"help", no_argument, NULL, 'h'},
                                      {"successor", no_argument, NULL, 's'},
                                      {NULL, 0, NULL, 0}};
  int opt_index = 0;
  while ((c = getopt_long(argc, argv, "hs", long_opts, &opt_index)) != -1) {
    switch (c) {
      case 'h':
        usage(argv[0]);
        return 0;
      case 's':
        args.successor = true;
        break;
      case '?':
        usage(argv[0]);
        return 1;
      default:
        printf("?? getopt returned character code 0%o ??\n", c);
    }
  }

  // Arg positional parsing
  if (argc - optind != 2) {
    std::cerr << "Error: Invalid number of arguments.\n\n";
    usage(argv[0]);
    return 1;
  } else {
    args.index_filename = argv[optind];
    args.kmers_filename = argv[optind + 1];
  }

  return 0;
};


int main(int argc, char* argv[]) {
  
  // Parse CLI args
  CliArgs args;
  args.successor = false;
  if (parse_args(argc, argv, args) != 0) {
    return 1;
  }

  // Deserialize index
  KmerManipulatorACGT km_index = KmerManipulatorACGT(31);
  auto index = ConwayBromageSD::deserialize(args.index_filename, &km_index);

  // Query index
  string line;
  std::ifstream kmers(args.kmers_filename);
  KmerManipulatorACGT km_query = KmerManipulatorACGT(30);
  if (args.successor) {

    int neighbors = 0;
    int total = 0;
    while (getline(kmers, line)) {
      auto query = km_query.encode(line);
      auto n = index.neighbours(query);
      if (n != 0) {
        neighbors++;
      }
      total++;
    }

    float_t percent = float(neighbors) / float(total);
    std::cout << neighbors << " / " << total << " (" << percent
              << "%) queried kmers had at least 1 neighbor." << std::endl;

  } else {
    int present = 0;
    int total = 0;

    while (getline(kmers, line)) {
      auto query = km_query.encode(line);
      present += index.contains(query);
      total++;
    }

    float_t percent = float(present) / float(total);
    std::cout << present << " / " << total << " (" << 100.0 * percent
              << "%) queried kmers present." << std::endl;
  }

  kmers.close();
}
