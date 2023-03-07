#include<bits/stdc++.h>
#include<unistd.h>
#include"tools.h"
using namespace std;
void decompmain(string& in_path,ofstream& fout)
{
    cout<<"decoding...\n";
    char ch;
    int ans_length,p=0;
    int cnt=0;
    ifstream fin1,fin2,fin3;
    fin1.open(in_path+".mine1");
    fin2.open(in_path+".mine2");
    fin3.open(in_path+".mine3");
    ans_length=in_int(fin1);
    rcnt=in_int(fin1);
    //cout<<ans_length<<'\n';
    
    while(p<ans_length){
        uint ch;
        ch=in_int(fin1);
        int cthre=p+16;//at most 4
        while(p<ans_length&&p<cthre){
            ans+=trans(int((ch&0xc0000000)>>30&3));//choose the highest pos 
            ch<<=2;
            p++;
        }
    }
    
    //cout<<ans<<'\n';
    
    vecr.resize(rcnt);
    int ppos=0,cpos,preread;
    for(int i=0;i<rcnt;i++){
        int bias,isrev,isrepeat,num;
        bias=in_char(fin2);
        cpos=bias+ppos;
        isrepeat=in_char(fin2);
        isrev=in_char(fin2);
        //cout<<"|||"<<i<<' '<<bias<<' '<<cpos<<' '<<isrepeat<<' '<<isrev<<'\n';
        if(isrepeat==1){
            if(isrev==0) vecr[i].str=vecr[preread].str;
            else vecr[i].str=cal_symm(vecr[preread].str);
        }
        else{
            preread=i;
            if(isrev==0){
                vecr[i].str=ans.substr(cpos,read_len);
            }
            else{
                for(int j=cpos+read_len-1;j>=cpos;j--){
                    vecr[i].str+=symm(ans[j],1);//ans has NO 'N',so use symm directly
                }
            }
            num=in_char(fin3);
            for(int j=0;j<num;j++){
                int dispos;
                char disc;
                dispos=in_char(fin3);
                disc=in_char(fin3);
                vecr[i].str[dispos]=disc;
            }
        }
        
        
        ppos=cpos;
    }
    
    if(order_preserve==0){
        for(int i=0;i<rcnt+repeatcnt;i++) fout<<vecr[i].str<<'\n';
    }
}