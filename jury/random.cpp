#include<bits/stdc++.h>
#include<unistd.h>
#include <sys/time.h>
using namespace std;
char trans(int c)//when compressing, use only 0-3
{
    char ret;
    switch(c){
        case 0:
            ret='A';
            break;
        case 1:
            ret='C';
            break;
        case 2:
            ret='G';
            break;
        case 3:
            ret='T';
            break;
        default:
            ret='N';
    }
    return ret;
}
signed main()
{
    srand(time(0));
    int n=100;
    int l=6;
    ofstream fout;
    fout.open("test.fastq");
    for(int i=1;i<=n;i++){
        fout<<i-1<<'\n';
        for(int j=0;j<l;j++){
            fout<<trans(rand()%5);
        }
        fout<<'\n';
        fout<<"x\n";
        fout<<"x\n";
    }
}