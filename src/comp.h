#include<bits/stdc++.h>
#include<unistd.h>
#include"tools.h"
#include"comp_prework.h"
#include"comp_align.h"
#include"comp_conmake.h"
#include"comp_realign.h"



void compmain(ifstream& fin)
{
    prework(fin);
    #if testflag
    for(int i=0;i<rcnt0;i++){
        cout<<i<<' '<<vread[i]<<' '<<isrepeat[i]<<' '<<repeatid[i]<<' '<<hasn[i]<<'\n';
    }
    #endif
    align();

    conmake();

    

    print_time();
    //test the length of ans
    // string tstr="";
    // for(int i=0;i<ccnt;i++){
    //     tstr+=vecc[i].str;
    // }
    // ofstream tfout;
    // tfout.open("test.tmp");
    // tfout<<tstr<<'\n';
    // cout<<tstr.length()<<'\n';
    // system("./bsc e test.tmp test.bsctmp");


    // //
    SCS_gen();
    #if testflag
    cout<<ans<<'\n';
    cout<<"info for all reads:\n";
    for(int i=0;i<rcnt0;i++){
        if(cpos[i]==-1){
            assert(isrepeat[i]);
            cout<<i<<' '<<"repeatid="<<repeatid[i]<<' '<<issymmrepeat[i]<<'\n';
        }
        else{
            cout<<i<<' '<<cpos[i]<<' '<<isrev[i]<<' '<<dismatch[i].size()<<'\n';
            for(auto it:dismatch[i]){
                cout<<it.first<<' '<<it.second<<'\n';
            }
        }
    }
    #endif
}