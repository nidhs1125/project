#include<bits/stdc++.h>
#include<unistd.h>
#include "tools.h"

int in_int(ifstream& fin)
{
    char ch;
    int ret=0;
    for(int i=0;i<4;i++){
        ret=ret<<8;
        fin.read(&ch,1);
        ret|=(unsigned char)ch;
    }
    return ret;
}

unsigned char in_char(ifstream& fin)
{
    char ch;
    fin.read(&ch,1);
    return (unsigned char)ch;
}


void decode(string& in_path,ofstream& fout)
{
    //system("./bsc d test.bsc1 test.mine1");
    //system("./bsc d test.bsc2 test.mine2");
    //system("./bsc d test.bsc3 test.mine3");
    cout<<"decoding...\n";
    char ch;
    int ans_length,p=0;
    int cnt=0;
    ifstream fin1,fin2,fin3;
    fin1.open(in_path+".mine1");
    fin2.open(in_path+".mine2");
    fin3.open(in_path+".mine3");
    ans_length=in_int(fin1);
    rcnt0=rcnt=in_int(fin1);
    read_len=in_int(fin1);
    cout<<ans_length<<' '<<rcnt<<' '<<read_len<<'\n';
    while(p<ans_length){
        ans+=in_char(fin1);
        p++;
    }
    #if testflag
    cout<<ans<<'\n';
    #endif
    vread.resize(rcnt);
    int ppos=0,cpos,preread;
    for(int i=0;i<rcnt;i++){
        string& str=vread[i];
        int bias,isrev,isrepeat,num;
        bias=in_char(fin2);
        cpos=bias+ppos;
        isrepeat=in_char(fin2);
        isrev=in_char(fin2);
        #if testflag
        cout<<"|||"<<i<<' '<<bias<<' '<<cpos<<' '<<isrepeat<<' '<<isrev<<'\n';
        #endif
        if(isrepeat==1){
            if(isrev==0) str=vread[preread];
            else str=cal_symm(vread[preread]);
        }
        else{
            preread=i;
            if(isrev==0){
                str=ans.substr(cpos,read_len);
            }
            else{
                string tmp=ans.substr(cpos,read_len);
                str=cal_symm(tmp);
            }
            num=in_char(fin3);
            for(int j=0;j<num;j++){
                int dispos;
                char disc;
                dispos=in_char(fin3);
                disc=in_char(fin3);
                str[dispos]=disc;
            }
        }
        #if testflag
        cout<<"+++"<<i<<' '<<str<<'\n';
        #endif
        
        ppos=cpos;
    }
    
    if(order_preserve==0){
        for(int i=0;i<rcnt+repeatcnt;i++) fout<<vread[i]<<'\n';
    }
}