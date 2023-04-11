#include<bits/stdc++.h>
#include<unistd.h>
#include"tools.h"

void realignment_prework(int k);
void realignment(int id,int pos,int is);//cur vecr id

void SCS_gen()
{
    cout<<"SCS_gen start\n";
    int cbias=0;
    for(int i=0;i<ccnt;i++){
        //if more than one, add
        ans+=vecc[i].str;
        vecc[i].spos=cbias;
        cbias+=vecc[i].str.length();
    }

    //for all read in ans , update the bias
    for(int i=0;i<rcnt0;i++){
        if(cid[i]!=-1){
            cpos[i]+=vecc[cid[i]].spos;
        }
    }
    //prework:hash all k-mer of reads not in ans

    #if testflag
    int k=4;
    #else
    int k=23;
    #endif

    realignment_prework(k);

    cout<<"realignment_prework finish\n";
    for(int i=0;i<(int)ans.length()-k;i++){//for all pos in ans
        ull tmp=get_seg_value(ans,i,i+k-1).first;
        if(mp_pos.find(tmp)==mp_pos.end()) continue;
        vector<srealign>& vec=mp_pos[tmp];
        for(srealign sr:vec){//for all id,check the dismatch,and update
            if(sr.isrev==0) realignment(sr.id,i-sr.pos,sr.isrev);
            else realignment(sr.id,i-(read_len-1-sr.pos),sr.isrev);
        }
    }

    //others
    int failednum=0;
    int prelen=ans.length();
    for(int i=0;i<rcnt0;i++){
        if(!isrepeat[i]&&cpos[i]==-1){
            failednum++;
            cpos[i]=ans.length();
            ans+=vread[i];
            dismatch[i].clear();
        }
    }
    int aftlen=ans.length();
    printf("failednum=%d,prelen=%d,aftlen=%d\n",failednum,prelen,aftlen);
    cout<<"SCS_gen finish\n";
}

void realignment_prework(int k)
{   
    assert(k<32);
    assert(read_len>=k);
    ull mask=(1ull<<2*k)-1;
    ull cur=0;

    int r1=max(read_len/2,k-1),l1=r1-k+1;
    int l2=min(read_len-k,read_len/2+1),r2=l2+k-1;
    int singletonnum=0;
    for(int i=0;i<rcnt0;i++){
        if(isrepeat[i]) continue;//repeat
        if(cpos[i]!=-1) continue;
        
        //else singleton
        puu tmp=get_seg_value(vread[i],l1,r1);
        mp_pos[tmp.first].pb(srealign{i,l1,0});
        mp_pos[tmp.second].pb(srealign{i,r1,1});
        tmp=get_seg_value(vread[i],l2,r2);
        mp_pos[tmp.first].pb(srealign{i,l2,0});
        mp_pos[tmp.second].pb(srealign{i,r2,1});
        singletonnum++;
    }
    printf("singletonnum=%d\n",singletonnum);
    
}
void realignment(int id,int pos,int is)
{
    if(pos<0||pos+read_len-1>ans.length()-1) return;
    int discnt=0;
    
    if(is==0){
        for(int i=0;i<read_len;i++){
            if(ans[pos+i]!=vread[id][i]) discnt++;
        }
    }
    else{
        for(int i=0;i<read_len;i++){
            if(symm(ans[pos+read_len-1-i],1)!=vread[id][i]) discnt++;
        }
    }
    if(discnt>maxdiscnt) return;
    if(cpos[id]==-1||dismatch[id].size()>discnt){
        dismatch[id].clear();
        isrev[id]=is;
        cpos[id]=pos;
        if(is==0){
            for(int i=0;i<read_len;i++){
                if(ans[pos+i]!=vread[id][i]) dismatch[id].pb(pic(i,vread[id][i]));
            }
        }
        else{
            for(int i=0;i<read_len;i++){
                if(symm(ans[pos+read_len-1-i],1)!=vread[id][i]) dismatch[id].pb(pic(i,vread[id][i]));
            }
        }
    }
}