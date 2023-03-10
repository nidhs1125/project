#include<bits/stdc++.h>
#include<unistd.h>
#include <sys/time.h>
#include"tools.h"
#include"comp.h"
#include"decomp.h"
using namespace std;

signed main(int argc,char** argv)
{
    timeval start, end; 
    gettimeofday(&start, NULL);
    
    char opt;  //
    const char *optstring = "i:o:dt:rp:";
    string in_path="",out_path="";
    bool flag=0;//flag:if compression;
    while ((opt = getopt(argc, argv, optstring)) != -1) {
        switch(opt){
            case 'i':
                in_path=optarg;
                break;
            case 'o':
                out_path=optarg;
                break;
            case 'd':
                flag=1;
                break;
            case 't':
                threshold=atoi(optarg);
                break;
            case 'r':
                order_preserve=1;
                break;
            case 'p':
                thread_num=atoi(optarg);
                break;
            default:
                cerr<<"Command line parameters default\n";
                exit(0);
        }
    }
    if(flag==0){
        if(in_path==""){
            #if testflag 
            in_path="test.fastq";
            #else
            in_path="../data/SRR554369_1.fastq";
            #endif
            //cerr<<"in_file lack\n";
            //exit(0);
        }
        if(out_path==""){
            out_path="test";
            //cerr<<"out_file lack\n";
            //exit(0);
        }
    }
    else{
        if(in_path==""){
            in_path="test";
        }
        if(out_path==""){
            out_path="decomp.fastq";
        }
    }
    
    
    
    if(flag==0){
        ifstream fin;
        fin.open(in_path,ios::in|ios::binary);
        if(!fin.is_open()){
            cerr<<"path default\n";
            exit(0);
        }
        compmain(fin,out_path);
    }
    else{
        ofstream fout;
        fout.open(out_path);
        if(!fout.is_open()){
            cerr<<"path default\n";
            exit(0);
        }
        decompmain(in_path,fout);
    }
    
    cout<<"finish\n";
    gettimeofday(&end, NULL);
    cout<<"time consumption:"<<1000*(end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)/1000<<"ms\n";
}
