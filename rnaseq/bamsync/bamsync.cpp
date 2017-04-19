//  bamsync: utility for patching a BAM file with FailedQC flags (0x200) from 
//  a reference BAM and adding read group IDs to each read based on read ID (RG:Z:xxxxx.x)
// 
//  Copyright (c) 2015 Francois Aguet, Broad Institute
// 
// 
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
// 
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
// 
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.


#include <iostream>
#include <iterator>
#include <vector>
#include <unordered_set>
#include <chrono>
#include <ctime>
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>

// Bamtools
#include "api/BamReader.h"
#include "api/BamWriter.h"

using namespace std;
using namespace BamTools;
namespace po = boost::program_options;


string timestamp() {
	time_t rawtime;
	char buffer[100];
	time(&rawtime);
	strftime(buffer, 100, "%H:%M:%S", localtime(&rawtime));
    string timestr = buffer;
    return "["+timestr+"]";
};


int main (int argc, char* argv[]) {
try {
        vector<string> infiles;
        string outfile;

        // option parser, help message
        po::options_description generic("Options");
        generic.add_options()
            ("help,h", "Print this message")
            ("version,v", "Display version information")
            ("output,o", po::value<string>(), "Specify output name. Default: <target_name>.patched.bam")
            ;

        po::options_description hidden("Hidden options");
        hidden.add_options()
            ("input-files", po::value< vector<string> >(), "Input files")
            ;

        po::options_description cmdline_options;
        cmdline_options.add(generic).add(hidden);

        po::positional_options_description p;
        p.add("input-files", -1);
   
        po::variables_map vm;
        store(po::command_line_parser(argc, argv).
              options(cmdline_options).positional(p).run(), vm);
        
        // check inputs
        if (vm.count("help")) {
            cout << "\n bamsync patches a target bam with FailedQC flags (0x200) from a reference bam\n and adds a read group ID to each read based on the read ID (RG:Z:xxxxx.x)\n\n" << generic << "\n" << "Usage:\n  bamsync reference.bam target.bam\n\n";
            return 0;
        }
   
        if (vm.count("version")) {
            cout << "bamsync v0.0.1\n";
            return 0;
        }
        
        if (!vm.count("input-files")) {
            cerr << "Inputs must be paths to 2 bam files. Run 'bamsync --help' for details.\n";
            return 0;
        } else {
            infiles = vm["input-files"].as< vector<string> >();
            if (infiles.size()!=2) {
                cerr << "Inputs must be paths to 2 bam files. Run 'bamsync --help' for details.\n";
                return 0;
            }
            
            if (vm.count("output")) {
                outfile = vm["output"].as< string >();
            } else {
                size_t d = infiles[1].rfind("/");
                size_t p = infiles[1].rfind(".");
                outfile = infiles[1].substr(d+1,p-d-1) + ".patched.bam";
            }
            
            BamReader reader1, reader2;
            if (!reader1.Open(infiles[0])) {
                cerr << "Could not open " << infiles[0] << endl;
                return 0;
            }
            if (!reader2.Open(infiles[1])) {
                cerr << "Could not open " << infiles[1] << endl;
                return 0;
            }
            
            // run
            cout << timestamp() << " Started sync." << endl;
            
            const SamHeader header1 = reader1.GetHeader();
            const SamHeader header2 = reader2.GetHeader();
            
            cout << timestamp() << " Parsing reference BAM for FailedQC reads (0x200 flag) ... " << flush;
            BamAlignment al;
            unordered_set<string> failedQC;
            while ( reader1.GetNextAlignment(al) ) {
                if (al.IsFailedQC()) {
                    failedQC.insert(al.Name);
                }
            }
            reader1.Close();
            cout << "done." << endl;
            
            // parse headers
            vector<string> hd1, rg1, pg1, co1;
            vector<string> hd2, sq2, rg2, pg2, co2;
                        
            boost::char_separator<char> sep("\n");
            string id;
            
            // source header lines
            string h1string = header1.ToString();
            boost::tokenizer< boost::char_separator<char> > tok(h1string, sep);
            for(auto it=tok.begin(); it!=tok.end(); ++it) {
                id = (*it).substr(0,3);
                if (id=="@HD") {
                    hd1.push_back(*it);
                // } else if (id=="@SQ") {
                //     sq1.push_back(*it);
                } else if (id=="@RG") {
                    rg1.push_back(*it);
                } else if (id=="@PG") {
                    pg1.push_back(*it);
                } else if (id=="@CO") {
                    co1.push_back(*it);
                }
            }
            
            // target header lines
            string h2string = header2.ToString();
            tok = boost::tokenizer< boost::char_separator<char> >(h2string, sep);
            for(auto it=tok.begin(); it!=tok.end(); ++it) {
                id = (*it).substr(0,3);
                if (id=="@HD") {
                    hd2.push_back(*it);
                } else if (id=="@SQ") {
                    sq2.push_back(*it);
                } else if (id=="@RG") {
                    rg2.push_back(*it);
                } else if (id=="@PG") {
                    pg2.push_back(*it);
                } else if (id=="@CO") {
                    co2.push_back(*it);
                }
            }

            bool is_genome_bam = true;
            if (hd2.size()==0) {  // transcriptome BAM: only has @SQ (for each transcript) and @RG lines
                is_genome_bam = false;
            }

            // new header:
            //   @HD line from target
            //   @SQ lines from target (may be different if reference changes)
            //   @RG lines from source
            //   @PG, @CO lines from target
            string header = "";
            if (is_genome_bam) {
                for (auto k=hd2.begin();k!=hd2.end();++k) {
                    header += *k + "\n";
                }
                for (auto k=sq2.begin();k!=sq2.end();++k) {
                    header += *k + "\n";
                }
                for (auto k=rg1.begin();k!=rg1.end();++k) {
                    header += *k + "\n";
                }
                for (auto k=pg2.begin();k!=pg2.end();++k) {
                    header += *k + "\n";
                }
                for (auto k=co2.begin();k!=co2.end();++k) {
                    header += *k + "\n";
                }
            } else {
                cout << "Number of transcript references (@SQ): " << sq2.size() << endl;
                for (auto k=sq2.begin();k!=sq2.end();++k) {
                    header += *k + "\n";
                }
                for (auto k=rg1.begin();k!=rg1.end();++k) {
                    header += *k + "\n";
                }
            }
            
            // chromosome names & lenghts
            const RefVector references = reader2.GetReferenceData();
                        
            cout << timestamp() << " Adding QC flag and read group ID to target BAM ... " << flush;
            BamWriter writer;
            if ( !writer.Open(outfile, header, references) ) {
                cerr << "Could not open output BAM file" << endl;
                return 0;
            }
            
            unordered_set<string>::const_iterator k;
            int i = 0;
            while ( reader2.GetNextAlignment(al) ) {
                k = failedQC.find(al.Name);
                // set QC flag
                if (k != failedQC.end()) {
                    al.SetIsFailedQC(true);
                    ++i;
                }
                
                // read group definition: RG:Z:${readID_comp_1[:5]}:${read_id_comp_2}
                size_t p = al.Name.find(":");
                al.EditTag("RG", "Z", al.Name.substr(0,5) + "." + al.Name.substr(p+1,1));
                writer.SaveAlignment(al);
            }                    
            writer.Close();
            reader2.Close();
            cout << "done." << endl;			
		    cout << timestamp() << " Finished sync. Number of QC-failed read mates (0x200 flag) patched: " << i << endl;
        }
    
    } catch(exception& e) {
        cout << "Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}
