#include<bits/stdc++.h>
#include<unistd.h>
#include "tools.h"



void encode(string& out_path)
{
    cout<<"encoding...\n";

    // for(auto it:vecr) cout<<it.rid<<' ';
    // cout<<'\n';
    //according to out_path, determine the out_path
    //has 3 paths in total: contig, readpos(bias) and the mismatch
    ofstream fout1,fout2,fout3;
    fout1.open(out_path+".mine1");
    fout2.open(out_path+".mine2");
    fout3.open(out_path+".mine3");

    //part 1
    out_int(fout1,ans.length());
    out_int(fout1,rcnt0);
    out_int(fout1,read_len);
    cout<<rcnt<<' '<<ccnt<<' '<<ans.length()<<'\n';
    //cout<<ans<<'\n';
    for(int i=0;i<ans.length();i++){
        out_char(fout1,ans[i]);
    }

    //part 2
    for(int i=0;i<rcnt0;i++){
        if(isrepeat[i]){
            cpos[i]=cpos[repeatid[i]];
            assert(cpos[i]!=-1);
        }
    }
    vecr.clear();
    for(int i=0;i<rcnt0;i++) vecr.pb(i);
    sort(vecr.begin(),vecr.end(),[](int id1,int id2){
        if(cpos[id1]!=cpos[id2]) return cpos[id1]<cpos[id2];
        else{
            int x,y;
            if(isrepeat[id1]==1) x=repeatid[id1];
            else x=id1;
            if(isrepeat[id2]==1) y=repeatid[id2];
            else y=id2;
            if(x!=y) return x<y;
            else return isrepeat[id1]<isrepeat[id2];
        }
    });
    int ppos=0;

    int dismatchcnt=0;
    int cntt=0;
    for(int j=0;j<rcnt0;j++){
        int i=vecr[j];
        out_char(fout2,cpos[i]-ppos);//bias
        // if(!(cpos[i]-ppos>=0&&cpos[i]-ppos<=read_len+1)){
        //     cout<<"***"<<i<<' '<<cpos[i]<<' '<<ppos<<'\n';
        // }
        //assert(cpos[i]-ppos>=0&&cpos[i]-ppos<=read_len+1);
        out_char(fout2,isrepeat[i]==1);//if repeat
        if(isrepeat[i]){
            out_char(fout2,issymmrepeat[i]);//if symm
        }
        else{
            out_char(fout2,(char)(isrev[i]));//if rev
            out_char(fout3,(char)dismatch[i].size());
            //may not satisfied,but few
            if(dismatch[i].size()>maxdiscnt){
                assert(pre[i]!=-1);
                cntt++;
            }
            dismatchcnt+=dismatch[i].size();
            if(dismatch[i].size()!=0){
                for(auto it1:dismatch[i]){
                    out_char(fout3,(char)it1.first);
                    out_char(fout3,(char)it1.second);
                }
            }
        }
        ppos=cpos[i];
    }
    cout<<dismatchcnt<<' '<<repeatcnt<<' '<<cntt<<'\n';

    fout1.close();
    fout2.close();
    fout3.close();
    string str="./bsc e "+out_path+".mine1 "+out_path+".bsc1";
    system(str.c_str());
    str="./bsc e "+out_path+".mine2 "+out_path+".bsc2";
    system(str.c_str());
    str="./bsc e "+out_path+".mine3 "+out_path+".bsc3";
    system(str.c_str());



    ofstream fout;
    ifstream fin;
    char ch;

    fout.open(out_path+".mine");//test.mine
    fin.open(out_path+".bsc1");
    while(fin>>ch) fout<<ch;
    fin.close();
    fin.open(out_path+".bsc2");
    while(fin>>ch) fout<<ch;
    fin.close();
    fin.open(out_path+".bsc3");
    while(fin>>ch) fout<<ch;
    fin.close();
    fout.close();
}