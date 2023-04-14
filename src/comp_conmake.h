#include<bits/stdc++.h>
#include<unistd.h>
#include"tools.h"


void* cal_block(void* arg);

//make contig
//the graph has built
void conmake()
{
    cout<<"contig_make start\n";


    vis.assign(rcnt0,0);
    cpos.assign(rcnt0,-1);
    dismatch.resize(rcnt0);
    isrev.assign(rcnt0,0);
    cid.assign(rcnt0,-1);

    int tmp;
    mutex_init();
    for(int i=0;i<0/*thread_num*/;i++){
        thread_tmp[i]=i;
        if((tmp=pthread_create(&thread_id[i],NULL,cal_block,&thread_tmp[i]))!=0){
            cout<<"pthread_create ERROR\n";
            exit(0);
        }
    }
    cal_block(&thread_tmp[thread_num]);
    for(int i=0;i<0/*thread_num*/;i++){
        void* tmp;
        pthread_join(thread_id[i],&tmp);
    }

    cout<<"contig make finish\n";
}

void* cal_block(void* arg)
{
    int id=*(int*)arg;
    while(1){
        int contig_id;
        int curid;
        pthread_mutex_lock(&my_mutex);
        while(my_id<rcnt0){
            if(bas_pre[my_id]==-1){
                my_id++;
                continue;
            }
            else if(find(bas_pre,my_id)==my_id) break;
            my_id++;
        }
        if(my_id>=rcnt0){
            pthread_mutex_unlock(&my_mutex);
            break;
        }
        curid=my_id++;
        if(block[curid].size()==1&&basket[block[curid][0]].size()==1){//only one read
            #if testflag
            printf("only one read:%d\n",my_id-1);
            #endif
            pthread_mutex_unlock(&my_mutex);
            continue;
        }
        else{
            contig_id=ccnt;
            vecc.pb(Contig());
            ccnt++;
            pthread_mutex_unlock(&my_mutex);
        }
        //for this block, make a contig
        #if testflag
        // cout<<"block_root="<<curid<<'\n';
        // for(int j=0;j<block[curid].size();j++){
        //     cout<<block[curid][j]<<' '<<block_rev[curid][j]<<' '<<block_bias[curid][j]<<' '<<vread[block[curid][j]]<<'\n';

        // }
        #endif



        int L=1e9,R=-1e9;
        string constr="";
        for(int i=0;i<block[curid].size();i++){
            L=min(L,block_bias[curid][i]);
            R=max(R,block_bias[curid][i]+read_len);
        }
        //printf("cal_block: %d %d %d\n",id,L,R);
        for(int j=L;j<R;j++){//for every pos
            int cnt[4]={0,0,0,0};
            for(int i=0;i<block[curid].size();i++){//for every basket
                for(int m=0;m<basket[block[curid][i]].size();m++){//for every read 
                    int id=basket[block[curid][i]][m];
                    int tmprev=basket_rev[block[curid][i]][m]^block_rev[curid][i];
                    int bias=block_bias[curid][i];
                    if(bias<=j&&read_len-1+bias>=j){
                        int curc;
                        if(tmprev==0) curc=trans(symm(vread[id][j-bias],0));
                        else curc=trans(symm(vread[id][read_len-1-(j-bias)],1));
                        //scout<<"debug: "<<j<<' '<<id<<' '<<tmprev<<' '<< bias<<' '<<curc<<'\n';
                        cnt[curc]++;
                    }
                }
            }
            //cout<<"++++"<<j<<' '<<cnt[0]<<' '<<cnt[1]<<' '<<cnt[2]<<' '<<cnt[3]<<'\n';
            int mx=cnt[3],mxpos=3;
            for(int i=2;i>=0;i--) if(cnt[i]>=mx) mx=cnt[i],mxpos=i;
            constr+=trans(mxpos);
        }
        //printf("|||cal_block: %d %d %d\n",id,L,R);
        for(int i=0;i<block[curid].size();i++){//for every basket
            for(int m=0;m<basket[block[curid][i]].size();m++){//for every read
                int id=basket[block[curid][i]][m];
                int tmprev=basket_rev[block[curid][i]][m]^block_rev[curid][i];
                int bias=block_bias[curid][i];
                cpos[id]=bias-L;
                isrev[id]=tmprev;
                cid[id]=contig_id;
                for(int j=0;j<vread[id].length();j++){
                    if(tmprev==0){
                        if(vread[id][j]!=constr[j+bias-L]){
                            dismatch[id].pb(pii(j,vread[id][j]));
                        }
                    }
                    else{//reverse
                        if(vread[id][j]!=symm(constr[read_len-1-j+bias-L],1)){//contig has no 'N'
                            dismatch[id].pb(pii(j,vread[id][j]));
                        }
                    }
                }
                if(dismatch[id].size()>maxdiscnt){//can't match the contig,the threshold can be more loose
                    //be singleton
                    #if testflag
                    printf("read be singleton:id=%d\n",id);
                    #endif
                    cpos[id]=-1;
                    isrev[id]=0;
                    cid[id]=-1;
                }
            }
        }
        pthread_mutex_lock(&array_mutex);
        vecc[contig_id].str=constr;
        pthread_mutex_unlock(&array_mutex);
    }
    return (void*)0;
}