
#include <ConwayBromageDatastructure.h>


void usage(std::string name) {
  std::cerr << "Build a Conway-Bromage structure from kmers and save it to disk.\n"
            << "Usage: " << name << " [OPTIONS] <kmers> <file_out>\n"
            << "\n"
            << "[ARGUMENTS]\n"
            << "    <kmers>    file containing sorted 31-mers to index\n"
            << "    <file_out> file to save the index to\n"
            << "\n"
            << "[OPTIONS]\n"
            << "    -h --help  Show this message" << std::endl;
};

int main(int argc, char* argv[]) {
  if (argc != 2 && argc != 3) {
    std::cerr << "Error: Invalid number of arguments. \n\n\n";
    usage(argv[0]);
    return 1;
  }

  std::string arg1;
  std::string arg2;

  if (argc == 2) {
    arg1 = argv[1];
    if (arg1 == "-h" || arg1 == "--help") {
      usage(argv[0]);
      return 0;
    } else {
      std::cerr << "Error: Invalid number of arguments. \n\n\n";
      usage(argv[0]);
      return 1;
    }
  }

  arg1 = argv[1];
  arg2 = argv[2];

  // Open kmer file
  std::ifstream file(arg1);

  // Build CBD index
  KmerManipulatorACGT km = KmerManipulatorACGT(31);
  ConwayBromageSD index(file, &km);
  file.close();

  // Save index to disk
  index.serialize(arg2);
}
