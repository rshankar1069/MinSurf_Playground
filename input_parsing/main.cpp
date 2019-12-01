#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include "atmsp.h"

int main(int argc, char* argv[])
{
    std::string filename,line;
    std::string Nx,Ny;
    std::string bottom,right,top,left;

    if(argc<1)
    {
        std::cout << "Please pass a file to be read with the program !!" << std::endl;
        abort();
    }
    else
        filename = argv[1];

    std::ifstream inp_file(filename);
    if(inp_file.is_open())
    {
        while(std::getline(inp_file,line))
        {
            if (line[0] == '#' || line == "" || line[0] == ' ')
                continue;
            else
            {
                std::vector<std::string> tokens;
                for(auto i = strtok(&line[0], " "); i != NULL; i = strtok(NULL, " "))
                    tokens.push_back(i);
                
                if(tokens[0] == "Nx")
                    Nx = std::stoi(tokens[1]);
                else if(tokens[0] == "Ny")
                    Ny = std::stoi(tokens[1]);
                else if(tokens[0] == "BC.BOTTOM")
                    bottom = tokens[1];
                else if(tokens[0] == "BC.RIGHT")
                    right = tokens[1];
                else if(tokens[0] == "BC.TOP")
                    top = tokens[1];
                else if(tokens[0] == "BC.LEFT")
                    left = tokens[1];
                else
                    break;
            }
        }
        inp_file.close();
    }
    else
    {
        std::cout << "Unable to open file" << std::endl;
    }

    /* Using the equation parser */

    ATMSP<double> parser;
    ATMSB<double> byteCode;
    // ATMSP<double> parserBottom,parserRight,parserTop,parserLeft;
    // ATMSB<double> byteCodeBottom,byteCodeRight,byteCodeTop,byteCodeLeft;

    std::vector<double> xval;

    int N = 100;
    double dt = 1.0/N;
    int i=0;

    while(i*dt <= 1.0)
    {
        xval.push_back(dt*i);
        i++;
    }

    parser.parse(byteCode,"sin(2*$pi*x)","x");
    i = 0;
    std::string writefile = "output.txt";
    std::ofstream out_file(writefile);
    while(i*dt <= 1.0)
    {
        byteCode.var[0] = xval[i];
        out_file << xval[i] << "\t" << byteCode.run() << std::endl;
        i++;
    }
    return 0;
}