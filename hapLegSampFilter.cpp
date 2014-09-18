#include <getopt.h>
#include <string>
#include <exception>
#include <boost/algorithm/string.hpp>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;
namespace io = boost::iostreams;

class site {
public:
  string id;
  unsigned pos;
  string a0;
  string a1;

  site(string i, unsigned p, string ia0, string ia1)
      : id(std::move(i)), pos(p), a0(std::move(ia0)), a1(std::move(ia1)) {};
};

void printFilterHapLeg(const string &keepSitesFile, const string &hapFile,
                       const string &legFile, const string &base) {

  // read in keep sites
  vector<unsigned> keepSites;
  clog << "Reading in sites to keep [" << keepSitesFile << "]" << endl;
  ifstream ksFS{ keepSitesFile };
  if (ksFS.is_open()) {
    string line;
    while (getline(ksFS, line)) {
      try {
        keepSites.push_back(stoul(line));
      }
      catch (const std::exception &e) {
        cerr << "Could not convert string to integer while reading keep sites "
                "file [" << keepSitesFile << "]: " << e.what() << endl;
      }
    }
  } else
    throw runtime_error("Could not open file " + keepSitesFile);

  // keep sites is no ready for binary lookup
  std::sort(keepSites.begin(), keepSites.end());

  vector<site> legend;
  clog << "Reading in legend file [" << legFile << "]" << endl;
  std::ifstream legFS(legFile, std::ios_base::binary);
  if (!legFS.is_open())
    throw runtime_error("Could not open file " + legFile);
  try {
    io::filtering_istream in;
    if (legFile.substr(legFile.find_last_of(".") + 1) == "gz")
      in.push(io::gzip_decompressor());
    in.push(legFS);

    vector<string> cols;
    string line;
    // throw away header
    getline(in, line);
    while (getline(in, line)) {
      boost::split(cols, line, boost::is_any_of(" "));
      try {
        legend.push_back(site(std::move(cols[0]), stoul(cols[1]),
                              std::move(cols[2]), std::move(cols[3])));
      }
      catch (const std::exception &e) {
        cerr << "Could not convert string to integer while reading legend "
                "file [" << legFile << "]: " << e.what() << endl
             << " Offending line was: " << line << endl;
      }
    }
  }
  catch (const io::gzip_error &e) {
    std::cout << e.what() << '\n';
  }

  // ok, now we can read in haplotypes and print them as we go along
  // assume haplotypes are gzipped
  clog << "Filtering haplotypes file [" << hapFile
       << "] and writing output files" << endl;
  ofstream outFS{ base + ".hap.gz", std::ios_base::binary };
  io::filtering_ostream out;
  out.push(io::gzip_compressor());
  out.push(outFS);

  ofstream outLeg{ base + ".legend" };
  if (!outLeg.is_open())
    throw runtime_error("could not open legend output file: " + base +
                        ".legend");
  else
    outLeg << "ID pos allele0 allele1\n";

  std::ifstream file{ hapFile, std::ios_base::binary };
  try {
    io::filtering_istream in;
    if (hapFile.substr(hapFile.find_last_of(".") + 1) == "gz")
      in.push(io::gzip_decompressor());
    in.push(file);

    unsigned lineNum{ 0 };
    unsigned keptSiteNum{ 0 };
    string line;
    while (getline(in, line)) {
      if (lineNum == legend.size())
        throw runtime_error("haps file contains more lines (" +
                            to_string(lineNum) + ") than the legend (" +
                            to_string(legend.size()) + ")");
      if (std::binary_search(keepSites.begin(), keepSites.end(),
                             legend[lineNum].pos)) {
        if (lineNum % 1000 == 0)
          clog << "\r" << keptSiteNum << "/" << lineNum;
        const site &lSite = legend[lineNum];
        out << line << "\n";
        outLeg << lSite.id << " " << lSite.pos << " " << lSite.a0 << " "
               << lSite.a1 << "\n";
        ++keptSiteNum;
      }
      ++lineNum;
    }
    if (lineNum != legend.size())
      throw runtime_error("legend contains more lines (" +
                          to_string(legend.size()) + " than haps file (" +
                          to_string(lineNum) + ")");
  }
  catch (const io::gzip_error &e) {
    std::cout << e.what() << '\n';
  }
}

int main(int argc, char **argv) {

  string hapFile, legFile, sampFile, keepPop, keepGroup, keepSex, keepSites;
  string base{ "out" };

  // Declare the supported options.
  po::options_description desc("Allowed options");
  desc.add_options()("help", "produce help message")("hap", po::value<string>(),
                                                     "impute2 haplotypes file")(
      "leg", po::value<string>(), "impute2 legend file")(
      "keepSites", po::value<string>(),
      "list of sites to keep")("base", po::value<string>(), "output file base");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
    cout << desc << "\n";
    return 1;
  }
  if (vm.count("hap"))
    hapFile = vm["hap"].as<string>();
  if (vm.count("leg"))
    legFile = vm["leg"].as<string>();
  if (vm.count("keepSites"))
    keepSites = vm["keepSites"].as<string>();
  if (vm.count("base"))
    base = vm["base"].as<string>();

  // make sure the correct input files are given
  if (keepSites.empty())
    throw runtime_error(
        "Only option right now is to keep sites. Please specify it.");
  if (!keepSites.empty())
    if (hapFile.empty() || legFile.empty())
      throw runtime_error("Need to specify hap and legend file");

  printFilterHapLeg(keepSites, hapFile, legFile, base);

  exit(0);
}
