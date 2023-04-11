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
    for(int i=0;i<thread_num;i++){
        thread_tmp[i]=i;
        if((tmp=pthread_create(&thread_id[i],NULL,cal_block,&thread_tmp[i]))!=0){
            cout<<"pthread_create ERROR\n";
            exit(0);
        }
    }
    cal_block(&thread_tmp[thread_num]);
    for(int i=0;i<thread_num;i++){
        void* tmp;
        pthread_join(thread_id[i],&tmp);
    }

    cout<<"contig make finish\n";
}

int dfs(vector<Edge>& vec,int rt,int cbias,int crev,int pa)
{
    assert(rt>=0&&rt<rcnt0);
    if(vis[rt]==1) return 0;
    int ret=basket[rt].size();
    vec.pb(Edge(rt,cbias,crev));
    vis[rt]=1;
    for(Edge& it:bas_edge[rt]){
        if(it.bias<=max_bias){
            assert(it.to==pa||vis[it.to]==0);//it is guaranteed that the graph is forest
            if(it.to!=pa) ret+=dfs(vec,it.to,cbias+it.bias,crev^it.rev,rt);
            
        }
    }
    return ret;
}

void* cal_block(void* arg)
{
    int id=*(int*)arg;
    while(1){
        int contig_id;
        vector<Edge> seq;
        pthread_mutex_lock(&my_mutex);
        while(my_id<rcnt0){
            if(isrepeat[my_id]||hasn[my_id]){
                my_id++;
                continue;
            }
            else if(find(my_id)==my_id&&vis[my_id]==0) break;
            my_id++;
        }
        if(my_id>=rcnt0){
            pthread_mutex_unlock(&my_mutex);
            break;
        }
        int tot=dfs(seq,my_id,0,first_rev[my_id],my_id);//the whole block should be taged
        my_id++;
        if(tot==1){//only one read
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
        int L=1e9,R=-1e9;
        string constr="";
        for(Edge it:seq){
            L=min(L,it.bias);
            R=max(R,it.bias+read_len);
        }
        //printf("cal_block: %d %d %d\n",id,L,R);
        for(int j=L;j<R;j++){//for every pos
            int cnt[4]={0,0,0,0};
            for(int i=0;i<seq.size();i++){//for every basket
                for(int m=0;m<basket[seq[i].to].size();m++){//for every read 
                    int id=basket[seq[i].to][m];
                    int tmprev=basket_rev[seq[i].to][m];
                    if(seq[i].bias<=j&&(int)vread[id].length()-1+seq[i].bias>=j){
                        int curc;
                        if(seq[i].rev^tmprev==0) curc=trans(symm(vread[id][j-seq[i].bias],0));
                        else curc=trans(symm(vread[id][read_len-1-(j-seq[i].bias)],1));
                        //cout<<"debug: "<<j<<' '<<i<<' '<<seq[i].bias<<' '<<(seq[i].rev^tmprev)<<' '<<curc<<'\n';
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
        for(int i=0;i<seq.size();i++){//for every basket
            for(int m=0;m<basket[seq[i].to].size();m++){//for every read
                int id=basket[seq[i].to][m];
                int tmprev=basket_rev[seq[i].to][m];
                cpos[id]=seq[i].bias-L;
                isrev[id]=tmprev^seq[i].rev;
                cid[id]=contig_id;
                for(int j=0;j<vread[id].length();j++){
                    if(seq[i].rev^tmprev==0){
                        if(vread[id][j]!=constr[j+seq[i].bias-L]){
                            dismatch[id].pb(pii(j,vread[id][j]));
                        }
                    }
                    else{//reverse
                        if(vread[id][j]!=symm(constr[read_len-1-j+seq[i].bias-L],1)){//contig has no 'N'
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
        #if testflag
        printf("cal_block: %d %s\n",contig_id,constr.c_str());
        #endif
    }
    return (void*)0;
}