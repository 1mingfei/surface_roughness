#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector> 

//read cfg file
typedef double real_t;
typedef char *  str_t;
typedef int    int_t;

/*
 struct atom_struct {
 int_t   size_N;
 str_t   name[size_N];
 real_t  coord[size_N][3];
 }
 */
int main(){
    real_t H0[3][3];
    real_t A_scale;
    real_t tmp_r0,tmp_r1,tmp_r2,tmp_r3;
    int_t  ii,jj;
    int_t  line_count,line_count1;
    int_t  N_atoms;
    int_t  entry_count;
    str_t  tmp_s;
    std::string line,word,word_1;

    std::ifstream myfile;
    myfile.open("example/small.cfg");

    getline(myfile,line);
    while (1)
    {
        if (sscanf(line.c_str(),"Number of particles = %d", &N_atoms))
        {
            std::cout<<N_atoms<<" atoms in the system."<<std::endl;
            break;
        }
        getline(myfile, line); //keep reading next line till end of the file
    }
    // clear error flags and reset seek position to 0
    myfile.clear();
    myfile.seekg(0);
    //declare mass and atomic position x
    real_t mass[N_atoms];
    real_t x[N_atoms][3];
    while (myfile.good())
    {
        getline(myfile, line); //keep reading next line till end of the file
        if (sscanf(line.c_str(),"A = %lf Angstrom (basic length-scale)", &A_scale))
        {
            std::cout<<"scale factor: "<< A_scale <<std::endl;
        }
        if (sscanf(line.c_str(),"H0(%d,%d) = %lf A", &ii, &jj, &tmp_r1))
        {
            H0[ii-1][jj-1]=tmp_r1;
        }
        if(sscanf(line.c_str(),"entry_count = %d",&ii))
        {
            entry_count=ii;
            std::cout<<entry_count<<" entries."<<std::endl;
            break;
        }
    }
    ii=0; jj=0;

    // clear error flags and reset seek position to 0
    myfile.clear();
    myfile.seekg(0);
    line_count=0;
    while (myfile.good())
    {
        getline(myfile, line);
        line_count+=1;
        std::cout<<line_count<<std::endl;
        std::istringstream iss(line);
        ii=0;
        jj=0;
        bool letter_flag=false;
        while(iss >> word) 
        {
            ii+=1;
            if (std::any_of(std::begin(word), std::end(word), ::isalpha))
            {
                letter_flag=true;
            }
        }
        if (ii==1 && !letter_flag)
        {
            sscanf(line.c_str(),"%lf",&tmp_r0);
            while (1)
            {
                getline(myfile, line);
                std::istringstream isss(line);
                int_t iii=0;
                bool letter_flag_1=false;
                while(isss >> word_1) 
                {
                    iii+=1;
                    if (std::any_of(std::begin(word_1), std::end(word_1), ::isalpha))
                    {
                        letter_flag_1=true;
                    }
                    std::cout<<iii<<std::endl;
                }
                line_count+=1;
                if (iii== 3)
                {
                    sscanf(line.c_str(),"%lf %lf %lf",&tmp_r1,&tmp_r2,&tmp_r3);
                    mass[jj]=tmp_r0;
                    x[jj][0]=tmp_r1;
                    x[jj][1]=tmp_r2;
                    x[jj][2]=tmp_r3;
                    std::cout<<"line"<<line_count<<std::endl;
                    std::cout<<jj<<"   "<<mass[jj]<<std::endl;
                    std::cout<<x[jj][0]<<" "<<x[jj][1]<<" "<<x[jj][2]<<std::endl;
                    jj+=1;
                }
                else if (iii==1 && !letter_flag_1)
                {
                    std::cout<<"break inner loop here!"<<std::endl;
                    break;
                }
            }
        }
        myfile.clear();
        myfile.seekg(line_count-1);
        line_count-=1;

    }
    /*
    std::cout<<"H matrix:"<<std::endl;
    for (int i=0;i<3;i++)
    {
        for (int j=0;j<3;j++)
        {
            std::cout<<H0[i][j]<<"   ";
        }
        std::cout<<std::endl;
    }
    std::cout<<"atom mass and positions:"<<std::endl;
    for  (int i=0;i<N_atoms;i++)
    {
        std::cout<<i<<"   "<<mass[i]<<std::endl;
        std::cout<<x[i][0]<<" "<<x[i][1]<<" "<<x[i][2]<<std::endl;
    }
    */
    return 0;
}
