
#include <ConwayBromageDatastructure.h>

void usage(std::string name) {
  std::cerr << "Loads a serialized Conway-Bromage structure from disk, "
               "for benchmarking purposes.\n"
            << "Usage: " << name << " [OPTIONS] <file_in>\n"
            << "\n"
            << "[ARGUMENTS]\n"
            << "    <file_in>  file containing the CBD\n"
            << "\n"
            << "[OPTIONS]\n"
            << "    -h --help  Show this message" << std::endl;
};


int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cerr << "Error: Invalid number of arguments. \n\n\n";
    usage(argv[0]);
    return 1;
  }

  std::string arg1;

  arg1 = argv[1];
  if (arg1 == "-h" || arg1 == "--help") {
    usage(argv[0]);
    return 0;
  }

  // Deserialize Index
  KmerManipulatorACGT km = KmerManipulatorACGT(31);
  auto index = ConwayBromageSD::deserialize(arg1, &km);
}
