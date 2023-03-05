#include<bits/stdc++.h>
#include<unistd.h>
#include"tools.h"
using namespace std;
void decompmain(ifstream& fin,ofstream& fout)
{
    cout<<"decoding...\n";
    char ch;
    int ans_length,p=0;
    int cnt=0;
    ans_length=in_int(fin);
    //cout<<ans_length<<'\n';
    
    while(p<ans_length){
        uint ch;
        ch=in_int(fin);
        int cthre=p+16;//at most 4
        while(p<ans_length&&p<cthre){
            ans+=trans(int((ch&0xc0000000)>>30&3));//choose the highest pos 
            ch<<=2;
            p++;
        }
    }
    
    //cout<<ans<<'\n';
    rcnt=in_int(fin);
    vecr.resize(rcnt);
    for(int i=0;i<rcnt;i++){
        int pos,len,isrev,num;
        pos=in_int(fin);
        len=in_char(fin);
        isrev=in_char(fin);
        //cout<<"|||"<<i<<' '<<pos<<' '<<len<<' '<<isrev<<'\n';
        if(isrev==0){
            vecr[i].str=ans.substr(pos,len);
        }
        else{
            for(int j=pos+len-1;j>=pos;j--){
                vecr[i].str+=symm(ans[j],1);//ans has NO 'N',so use symm directly
            }
        }
        num=in_char(fin);
        for(int j=0;j<num;j++){
            int dispos;
            char disc;
            dispos=in_char(fin);
            disc=in_char(fin);
            vecr[i].str[dispos]=disc;
        }
    }
    repeatcnt=in_int(fin);
    vecr.resize(rcnt+repeatcnt);
    for(int i=rcnt;i<rcnt+repeatcnt;i++){
        int pos,isrev;
        pos=in_int(fin);
        isrev=in_char(fin);
        if(isrev==0) vecr[i].str=vecr[pos].str;
        else vecr[i].str=cal_symm(vecr[pos].str);
        if(order_preserve==1){
            vecr[i].rid=in_int(fin);
        }
    }
    if(order_preserve==0){
        for(int i=0;i<rcnt+repeatcnt;i++) fout<<vecr[i].str<<'\n';
    }
    else{
        int p=0,q=rcnt;
        int curid=0;
        while(curid!=rcnt+repeatcnt){
            if(q==rcnt+repeatcnt) vecr[p++].rid=curid++;
            else if(vecr[q].rid==curid){
                q++;
                curid++;
            }
            else vecr[p++].rid=curid++;
        }
        sort(vecr.begin(),vecr.end(),cmp4);
        for(int i=0;i<rcnt+repeatcnt;i++){
            fout<<vecr[i].str<<'\n';
        }
    }
}