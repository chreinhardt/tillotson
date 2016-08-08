#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

extern "C" {
#include "tipsy.h"
}

using std::cout;
using std::endl;
using std::cin;

struct particle
{
    double R;
    double rho;
    double T;
    double M;
};

class profile
{
private:
    double pi;
    double G;

public:
    double R;
    double Mtot;
    std::string datafilename;

    std::vector<particle> data;

    profile(std::string File)
    {
        datafilename = File;

        Mtot = 0;
        R = 0;

        pi = 3.14;
        G = 1;

    }

    void read_datafile()
    {
        std::ifstream datafile;

        //datafile.exceptions(std::ifstream::failbit | std::ifstream::badbit);

        particle p;

        datafile.open(datafilename.c_str());


        while(!datafile.eof())
        {
            datafile >> p.R;
            datafile >> p.rho;
            datafile >> p.T;
            datafile >> p.M;

            data.push_back(p);
        }
        datafile.close();
        Mtot = p.M;
        R = p.R;

    }


    double get_rho(double r)
    {
        int a,b;

        int c;

        a = 0;
        b = data.size()-1;

        while((b - a) > 1)
        {
            c = (a + b)/2;

            if(data[c].R >= r)
            {
                b = c;
            }
            if(data[c].R <= r)
            {
                a = c;
            }
        }

        double rho;

        rho = (data[b].rho + data[a].rho)/2.0;

        return rho;
    }

    double get_R(double M)
    {

        int a,b;

        int c;

        a = 0;
        b = data.size()-1;

        //cout << "a: " << a << " b: " << b << endl;

        while((b - a) > 1)
        {
            c = (a + b)/2;

            if(data[c].M >= M)
            {
                b = c;
            }
            if(data[c].M <= M)
            {
                a = c;
            }
        }

        double R;

        R = (data[b].R + data[a].R)/2.0;

        return R;
    }
};

int main(int argc, char** argv)
{
        unsigned int N;

        double dM, T;

        double x,R,rho,theta,phi;
        TCTX outTipsy = NULL;
        struct gas_particle gp;

        std::string Profile;
        std::string File;

        if (argc != 4) {
            fprintf(stderr,"Usage: writeprofile <profile> <N> <tipsyoutputfile.std>\n");
            exit(1);
        }

        Profile = argv[1];
        File = argv[2];

        N = atoi(argv[3]);
        assert(N > 0 && "N has to be larger than zero");

        profile p(Profile);

        cout << "Loading profile " << Profile << endl;
/*
        try {
            p.read_datafile();
        }
        catch (const std::ifstream::failure& e)
        {
            cout << "Error while opening/reading file: " << e.what() << endl;
            return 1;
        }
*/
        p.read_datafile();

        cout << "Profile loaded ..." << endl;

        cout << "M: " << p.Mtot << " R: " << p.R << endl;

        dM = p.Mtot / N;

        cout << "Creating initial conditions: " << File << " N: " << N << " dM: " << dM << endl;

        TipsyInitialize(&outTipsy,0,NULL);
        for (unsigned int i(0);i<N;++i) {
                R = p.get_R((i+1)*dM);

                x = 2.0*rand()/(RAND_MAX+1.0)-1.0;
                theta = acos(x);

                phi = 2.0*M_PI*rand()/(RAND_MAX+1.0);

                rho = p.get_rho(R);
                T = 0.001;

                cout << "R: " << R << " rho: " << rho << " M: " << (i+1)*dM << endl;
                gp.mass = dM;
                gp.pos[0] = R*sin(theta)*cos(phi);
                gp.pos[1] = R*sin(theta)*sin(phi);
                gp.pos[2] = R*cos(theta);
                gp.vel[0] = 0;
                gp.vel[1] = 0;
                gp.vel[2] = 0;
                gp.rho = rho;
                gp.temp = T;
                gp.hsmooth = 0.1;
                gp.metals = 0;
                gp.phi = 0;
                TipsyAddGas(outTipsy,&gp);
        }
        TipsyWriteAll(outTipsy,0.0,const_cast<char*>(File.c_str()));
        return 0;
}
