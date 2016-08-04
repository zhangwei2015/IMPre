#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<ctype.h>
#include<zlib.h>
#include "myhash.h"
#define MAX_FILENAME_LEN 256
#define INCREMENT 16
#define MAX_LEN 32
typedef unsigned char bit8_t;
bit8_t nst_nt4_table[256] = {
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};
char InFileName[MAX_FILENAME_LEN]="";
char OutDir[MAX_FILENAME_LEN]="";
char PrefixName[MAX_FILENAME_LEN]="";
char TempFileName[MAX_FILENAME_LEN]="";
char GeneType[10]="";
int SeedLen,ISeedLen,MaxSeedNum,MinSeedNum,FMaxSeedNum,MinSeqNum,VMaskNum,JMaskNum;
struct seq_t{
    int id;
	int num;
    char *seq;
};
struct seq_t *SeqInfo=NULL;
typedef struct{
	int len;
    bit64_t *val;
}nst_bfa1_t;
struct Seed_Value{
    char *seq;
	int len;
    bit64_t *val;
};
typedef struct{
    int len;
    bit64_t val[];
}lookup_key_t;
typedef struct SEED_TYPE{
	UT_hash_handle hh;
	int num;
	int len;
	bit64_t val[];
}seed_type;
seed_type *SeedInfo=NULL,*SeedOther=NULL,*SeedInfo_New=NULL;;//Seed information
struct ID_type{
	int seq_id;
	int pos;
	int read_num;
};
typedef struct BACK_TYPE{
		UT_hash_handle hh;
		int num,maxnum;
		struct ID_type *IDInfo;
		int len;
		bit64_t val[];
}back_type;
back_type *BackInfo=NULL,*BackInfo_New=NULL;//sequence id information

void Usage(void) {
    fprintf(stdout, "\nUsage:\tSeed Cluster software [options]\n");
    fprintf(stdout, "\t-i  <str>   Input file,*.fq or *.fq.gz\n");
    fprintf(stdout, "\t-o  <str>   Output directory\n");
    fprintf(stdout, "\t-g  <char>  Gene type, V or J\n");
    fprintf(stdout, "\t-l  <int>   Seed length[40]\n");
    fprintf(stdout, "\t-n  <int>   The maximum number of seed for output [400]\n");
	fprintf(stdout, "\t-s  <int>   The minimum supporting sequence for seed [2]\n");
	fprintf(stdout, "\t-x  <int>   The number of seeds considered for first spliting the sequences [-n*10]\n");
	fprintf(stdout, "\t-d  <int>   The minimum of unique sequence number for supporting a seed [2]\n");
	fprintf(stdout, "\t-p  <str>   The prefix name for output file\n");
	fprintf(stdout, "\t-V  <int>   The number of 3'-terminal bases masked for V [10]\n");
	fprintf(stdout, "\t-J  <int>   The number of 5'-terminal bases masked for J [5]\n");
    exit(EXIT_FAILURE);
}
/*Initialize Parameters*/
void optInit(int argc,char **argv)
{
    char c;
    while((c=getopt(argc,argv,"i:o:g:l:n:s:x:d:p:V:J:h"))!=-1)
        switch(c)
        {
            //essential arguments
            case 'i':snprintf(InFileName,MAX_FILENAME_LEN,"%s",optarg);break;
            case 'o':snprintf(OutDir,MAX_FILENAME_LEN,"%s",optarg);break;
            case 'g':snprintf(GeneType,10,"%s",optarg);break;
            case 'l':SeedLen=atoi(optarg);break;
            case 'n':MaxSeedNum=atoi(optarg);break;
			case 's':MinSeedNum=atoi(optarg);break;
			case 'x':FMaxSeedNum=atoi(optarg);break;
			case 'd':MinSeqNum=atoi(optarg);break;
			case 'p':snprintf(PrefixName,MAX_FILENAME_LEN,"%s",optarg);break;
			case 'V':VMaskNum=atoi(optarg);break;
			case 'J':JMaskNum=atoi(optarg);break;
            case 'h':
            case '?':
                Usage();
        }
}
/*Delete Return char in line*/
void chomp(const char *s)
{
    char *p;
    while (NULL != s && NULL != (p = strrchr(s,'\n'))){
        *p = '\0';
    }
}
nst_bfa1_t *nst_new_bfa1()
{
    nst_bfa1_t *bfa1;
    bfa1 = (nst_bfa1_t*)calloc(1,sizeof(nst_bfa1_t));
    bfa1->len = 0;
    return bfa1;
}
void nst_delete_bfa1(nst_bfa1_t *bfa1)
{
    if (bfa1 == 0) return;
    free(bfa1->val);
    free(bfa1);
}
nst_bfa1_t *StringToULL(char *str)
{
    int i,len,new_len;
    bit64_t s;
    nst_bfa1_t *bfa1;

    len=strlen(str);
    if(len<0) return NULL;
    bfa1 = nst_new_bfa1();
    bfa1->len = len>>5;
    if (len&0x1f) ++(bfa1->len);
    bfa1->val = (bit64_t*)malloc(sizeof(bit64_t) * bfa1->len);
    s=0ull;
    for (i = 0; i != len; ++i) {
        int tmp = nst_nt4_table[(int)str[i]];
        s <<= 2;
        if (tmp < 4) s |= tmp;
        if((i&0x1f) ==0x1f)bfa1->val[i>>5] = s;
    }
    if (len&0x1f){
        s <<= (32 - (i&0x1f)) << 1;
        bfa1->val[len>>5] = s;
    }
    return bfa1;
}
int compare1(const void *a,const void *b)
{
    struct seq_t *ia = (struct seq_t *)a;
    struct seq_t *ib = (struct seq_t *)b;
    if(ia->id>ib->id)return 1;
    if(ia->id==ib->id) return 0;
    if(ia->id<ib->id)return -1;
}
int compare2(const void *a,const void *b)
{
    seed_type *ia = (seed_type *)a;
    seed_type *ib = (seed_type *)b;
    if(ia->num>ib->num)return -1;
    if(ia->num==ib->num) return 0;
    if(ia->num<ib->num)return 1;
}
int compare3(const void *p, const void *q)
{
    return (*(int *)p - *(int *)q);
}
int compare4(const void *a,const void *b)
{
    int i,mark;
    struct Seed_Value *ia=(struct Seed_Value*)a;
    struct Seed_Value *ib=(struct Seed_Value*)b;
    if(ia->len>ib->len) return 1;
    if(ia->len==ib->len){
        mark=0;
        for(i=0;i<ia->len;i++){
            if(ia->val[i]==ib->val[i])continue;
            if(ia->val[i]>ib->val[i]){mark=1;break;}
            if(ia->val[i]<ib->val[i]){mark=-1;break;}
        }
        return mark;
    }
    if(ia->len<ib->len) return -1;
}
/*Fetch sequence number in input file*/
int Get_Seq_Num(char *id)
{
	char *id_t=NULL,*str_tok=NULL;
	int cyc_num,num;
	
	id_t=(char*)calloc(strlen(id)+1,sizeof(char));
	strcpy(id_t,id);
	str_tok=strtok(id_t,":");cyc_num=0;
   	while(str_tok!=NULL){
		if(cyc_num==1)num=atoi(str_tok);
		str_tok=strtok(NULL,":");cyc_num++;
    }
	if(id_t){free(id_t);id_t=NULL;}
	return num;
}
/*Fetch sequence id in input file*/
int Get_Seq_ID(char *id)
{
	char *id_t=NULL,*str_tok=NULL,*id_tmp=NULL;
    int cyc_num,num;

    id_t=(char*)calloc(strlen(id)+1,sizeof(char));
    strcpy(id_t,id);
    str_tok=strtok(id_t,":");cyc_num=0;
    while(str_tok!=NULL){
        if(cyc_num==0){
			id_tmp=(char*)calloc(strlen(str_tok),sizeof(char));
			strcpy(id_tmp,str_tok+2);
			num=atoi(id_tmp);
			if(id_tmp){free(id_tmp);id_tmp=NULL;}
			break;
		}
        str_tok=strtok(NULL,":");cyc_num++;
    }
    if(id_t){free(id_t);id_t=NULL;}
    return num;
}
/*Output Seed information*/
void Output_Seed(int Readnum,int Uniq_Seed_Num)
{
    FILE *fp;
    char *OutFileName=NULL,FinalSeedNum[10],command[MAX_FILENAME_LEN];
    struct Seed_Value *Seed_V,*Seed_P,*Seed_T;
    int i,j,index=0;
    struct seq_t *Find_P;
    seed_type *seed_tmp;
    back_type *back_tmp;
		

    /*Read temp value and seed corresponding file*/
    if((fp=fopen(TempFileName,"r"))==NULL){
        fprintf(stderr,"Can't open the seed and its value corresponding file!\n");
        exit(EXIT_FAILURE);
    }
    Seed_V=(struct Seed_Value*)calloc(Uniq_Seed_Num,sizeof(struct Seed_Value));
    for(i=0;i<Uniq_Seed_Num;i++){
        Seed_V[i].seq=(char*)calloc(SeedLen+1,sizeof(char));
        fread(Seed_V[i].seq,sizeof(char),SeedLen,fp);
        fread(&Seed_V[i].len,sizeof(int),1,fp);
        Seed_V[i].val=(bit64_t*)calloc(Seed_V[i].len,sizeof(bit64_t));
        for(j=0;j<Seed_V[i].len;j++)fread(&Seed_V[i].val[j],sizeof(bit64_t),1,fp);
    }
    fclose(fp);
    /*Remove temp file*/
    sprintf(command,"rm -f %s",TempFileName);
    system(command);
	/*output seed files*/
	unsigned keylen;
    qsort(Seed_V,Uniq_Seed_Num,sizeof(struct Seed_Value),compare4);
    HASH_SORT(SeedInfo_New,compare2);
    for(seed_tmp=SeedInfo_New;seed_tmp!=NULL;seed_tmp=(seed_type*)(seed_tmp->hh.next)){
        if(seed_tmp->num==0)continue;
		keylen=seed_tmp->len*sizeof(bit64_t);
        HASH_FIND(hh,BackInfo_New,seed_tmp->val,keylen,back_tmp);
        if(back_tmp){
            index++;
            OutFileName=(char*)calloc(strlen(OutDir)+50,sizeof(char));
            strcpy(OutFileName,OutDir);strcat(OutFileName,"/seed.");
            if(strlen(PrefixName)>0){strcat(OutFileName,PrefixName);strcat(OutFileName,".");}
            sprintf(FinalSeedNum,"%d",index);strcat(OutFileName,FinalSeedNum);
            strcat(OutFileName,".seq");
            if((fp=fopen(OutFileName,"w"))==NULL){
                fprintf(stderr,"Can't open the output file!\n");
                exit(EXIT_FAILURE);
            }
			Seed_T=(struct Seed_Value*)malloc(sizeof(sizeof(struct Seed_Value)));
			Seed_T->len=seed_tmp->len;Seed_T->val=(bit64_t*)calloc(seed_tmp->len,sizeof(bit64_t));
			for(i=0;i<seed_tmp->len;i++)Seed_T->val[i]=seed_tmp->val[i];
            Seed_P=(struct Seed_Value*)bsearch(Seed_T,Seed_V,Uniq_Seed_Num,sizeof(struct Seed_Value),compare4);
			free(Seed_T->val);Seed_T->val=NULL;
			free(Seed_T);Seed_T=NULL;
            if(Seed_P)fprintf(fp,"%s:%d:%d\n",Seed_P->seq,seed_tmp->num,back_tmp->num);
            for(i=0;i<back_tmp->num;i++){
                Find_P=(struct seq_t*)bsearch(&back_tmp->IDInfo[i].seq_id,SeqInfo,Readnum,sizeof(struct seq_t),compare1);
                if(Find_P)fprintf(fp,">%c%d:%d-%d:%d\n%s\n",tolower(GeneType[0]),Find_P->id,Find_P->num,back_tmp->IDInfo[i].pos,back_tmp->IDInfo[i].pos+SeedLen-1,Find_P->seq);
            }
            fclose(fp);
            sprintf(command,"gzip -f %s",OutFileName);
            system(command);
            if(OutFileName){free(OutFileName);OutFileName=NULL;}
        }
    }
    if(Seed_V){
		for(i=0;i<Uniq_Seed_Num;i++){free(Seed_V[i].val);Seed_V[i].val=NULL;free(Seed_V[i].seq);Seed_V[i].seq=NULL;}
		free(Seed_V);Seed_V=NULL;
	}
}	
/*Cluster Seed from each read*/
void Cluster_Seed(int Readnum)
{
	/*Define parameters*/
	int i,j,k,l;
	int min,max,seq_len;//Effective start and end positions for fetch seed 
	char *seed;//Seed sequence
	nst_bfa1_t *seed_v;//Seed value
	lookup_key_t *lookup_key;
	struct ID_type *Temp_P=NULL;
	FILE *fp;
	unsigned keylen;
	seed_type *seed_tmp;
	back_type *back_tmp;

	strcpy(TempFileName,OutDir);strcat(TempFileName,"/Seed_Value.txt");
	if((fp=fopen(TempFileName,"w"))==NULL){
		fprintf(stderr,"Can't open the seed and its value corresponding file!\n");
        exit(EXIT_FAILURE);		
	}
	int Uniq_Seed_Num=0;
	for(i=0;i<Readnum;i++){
		seq_len=strlen(SeqInfo[i].seq);
		min=0;max=seq_len-SeedLen;
		if(!strcmp(GeneType,"V"))max=max-VMaskNum;
		if(!strcmp(GeneType,"J"))min=JMaskNum;
		if(seq_len<SeedLen)continue;
		for(j=0;j<=seq_len-SeedLen;j++){
			seed=(char*)calloc(SeedLen+1,sizeof(char));
            strncpy(seed,SeqInfo[i].seq+j,SeedLen);
			seed_v=StringToULL(seed);
			lookup_key=(lookup_key_t*)malloc(sizeof(lookup_key_t)+seed_v->len*sizeof(bit64_t));
			memset(lookup_key,0,sizeof(lookup_key_t)+seed_v->len*sizeof(bit64_t));
			lookup_key->len = seed_v->len;
    		memcpy(lookup_key->val,seed_v->val,seed_v->len*sizeof(bit64_t));
			nst_delete_bfa1(seed_v);
			keylen=lookup_key->len*sizeof(bit64_t);
			if(j>=min && j<=max){
				/*Store seed information*/
				HASH_FIND(hh,SeedInfo,lookup_key->val,keylen,seed_tmp);
                if(!seed_tmp){
					seed_tmp=(seed_type*)malloc(sizeof(seed_type)+lookup_key->len*sizeof(bit64_t));
					memset(seed_tmp,0, sizeof(seed_type)+lookup_key->len*sizeof(bit64_t));
                    seed_tmp->len=lookup_key->len;
					memcpy(seed_tmp->val,lookup_key->val,lookup_key->len*sizeof(bit64_t));
                    seed_tmp->num=SeqInfo[i].num;
					HASH_ADD(hh,SeedInfo,val,keylen,seed_tmp);
					seed_tmp=NULL;
                }else seed_tmp->num+=SeqInfo[i].num;
			}else{
				/*Store seed information in mask region*/
                HASH_FIND(hh,SeedOther,lookup_key->val,keylen,seed_tmp);
                if(!seed_tmp){
					seed_tmp=(seed_type*)malloc(sizeof(seed_type)+lookup_key->len*sizeof(bit64_t));
                    memset(seed_tmp,0, sizeof(seed_type)+lookup_key->len*sizeof(bit64_t));
                    seed_tmp->len=lookup_key->len;
                    memcpy(seed_tmp->val,lookup_key->val,lookup_key->len*sizeof(bit64_t));
                    seed_tmp->num=SeqInfo[i].num;
                    HASH_ADD(hh,SeedOther,val,keylen,seed_tmp);
                    seed_tmp=NULL;
                }else seed_tmp->num+=SeqInfo[i].num;
			}
			/*Store ID information*/
			HASH_FIND( hh,BackInfo,lookup_key->val,keylen,back_tmp);
			if(!back_tmp){
				fwrite(seed,sizeof(char),SeedLen,fp);fwrite(&lookup_key->len,sizeof(int),1,fp);
				for(k=0;k<lookup_key->len;k++)fwrite(&lookup_key->val[k],sizeof(bit64_t),1,fp);
				Uniq_Seed_Num++;
				back_tmp=(back_type*)malloc( sizeof(back_type)+lookup_key->len*sizeof(bit64_t));
                memset(back_tmp,0, sizeof(back_type)+lookup_key->len*sizeof(bit64_t));
                back_tmp->len=lookup_key->len;
                memcpy(back_tmp->val,lookup_key->val,lookup_key->len*sizeof(bit64_t));
				back_tmp->num=back_tmp->maxnum=0;
				back_tmp->IDInfo=(struct ID_type*)malloc(sizeof(struct ID_type));//memset(back_tmp->IDInfo,0,sizeof(struct ID_type));
				back_tmp->IDInfo[back_tmp->num].seq_id=SeqInfo[i].id;
				back_tmp->IDInfo[back_tmp->num].pos=j+1;
				back_tmp->IDInfo[back_tmp->num].read_num=SeqInfo[i].num;back_tmp->num++;
				HASH_ADD(hh,BackInfo,val,keylen,back_tmp);
				back_tmp=NULL;
			}else{
				if(back_tmp->num+1>=back_tmp->maxnum){
					back_tmp->maxnum+=INCREMENT; 
					Temp_P=(struct ID_type*)realloc(back_tmp->IDInfo,back_tmp->maxnum*sizeof(struct ID_type));
					if(!Temp_P){
						fprintf(stderr,"Unable to allocate memory\n");
				        exit(EXIT_FAILURE);				
					}else back_tmp->IDInfo=Temp_P;
				}
				back_tmp->IDInfo[back_tmp->num].seq_id=SeqInfo[i].id;
				back_tmp->IDInfo[back_tmp->num].pos=j+1;
				back_tmp->IDInfo[back_tmp->num].read_num=SeqInfo[i].num;
				back_tmp->num++;
			}
			free(seed);seed=NULL;free(lookup_key);lookup_key=NULL;
		}/*End j & one sequence*/
	}/*End i & all sequence*/
	fclose(fp);
	printf("End Read Seed, The number of unique seed %d\n",HASH_COUNT(SeedInfo));
	printf("End Read Seed, The number of unique seed %d\n",HASH_COUNT(SeedOther));
	printf("End Read Seed, The number of unique seed %d\n",HASH_COUNT(BackInfo));
	/*for(seed_tmp=SeedInfo;seed_tmp!=NULL;seed_tmp=(seed_type*)(seed_tmp->hh.next)){
		printf("%d:\t",seed_tmp->key.len);
		for(i=0;i<seed_tmp->key.len;i++)printf("%llu\t",seed_tmp->key.val[i]);
		printf("%d\n",seed_tmp->num);
	}*/
	/*Check Seed in Mask region*/
    seed_type *seed_temp;
    for(seed_tmp=SeedOther;seed_tmp!=NULL;seed_tmp=(seed_type*)(seed_tmp->hh.next)){
		keylen=seed_tmp->len*sizeof(bit64_t);
        HASH_FIND(hh,SeedInfo,seed_tmp->val,keylen,seed_temp);
        if(seed_temp)seed_temp->num+=seed_tmp->num;
    }
	/*Free SeedOther*/
	seed_type *temp;
	HASH_ITER(hh,SeedOther,seed_tmp,temp) {
      HASH_DEL(SeedOther,seed_tmp);
      free(seed_tmp);
    }
	/*Delete ineffective id information*/
    for(back_tmp=BackInfo;back_tmp!=NULL;back_tmp=(back_type*)(back_tmp->hh.next)){
		keylen=back_tmp->len*sizeof(bit64_t);
        HASH_FIND(hh,SeedInfo,back_tmp->val,keylen,seed_temp);
        if(!seed_temp){HASH_DEL(BackInfo,back_tmp);free(back_tmp);}
    }
    /*Fetch the certain number of seeds based on descending order of reads number supporting seed*/
    HASH_SORT(SeedInfo,compare2);
    int Stat_Seed_Num=1;
    for(seed_tmp=SeedInfo;seed_tmp!=NULL;seed_tmp=(seed_type*)(seed_tmp->hh.next)){
        if(Stat_Seed_Num>FMaxSeedNum || seed_tmp->num<MinSeedNum){
			keylen=seed_tmp->len*sizeof(bit64_t);
            HASH_FIND(hh,BackInfo,seed_tmp->val,keylen,back_tmp);
            if(back_tmp){HASH_DEL(BackInfo,back_tmp);free(back_tmp);}
            HASH_DEL(SeedInfo,seed_tmp);free(seed_tmp);
        }else Stat_Seed_Num++;
    }
	printf("Starting filtering\n");
	/*Filter seed information*/
    int flag=0,index,New_Num;
    int *ID_arr=NULL,*Find_P,*Flag_ID;
    back_type *back_ttmp;
    int HASH_SIZE;
	unsigned keylen1,keylen2;

    HASH_SIZE=HASH_COUNT(SeedInfo);
    while(HASH_SIZE>0 && flag<MaxSeedNum){
        HASH_SORT(SeedInfo,compare2);
        seed_tmp=SeedInfo;
		keylen1=seed_tmp->len*sizeof(bit64_t);
        HASH_FIND(hh,BackInfo,seed_tmp->val,keylen1,back_tmp);
		if(back_tmp){
            if(seed_tmp->num>=MinSeedNum && back_tmp->num>=MinSeqNum){
				keylen2=back_tmp->len*sizeof(bit64_t);
                HASH_FIND(hh,BackInfo_New,back_tmp->val,keylen2,back_ttmp);
                if(!back_ttmp){
					 back_ttmp=(back_type*)malloc( sizeof(back_type)+back_tmp->len*sizeof(bit64_t));
                	memset(back_ttmp,0, sizeof(back_type)+back_tmp->len*sizeof(bit64_t));
                	back_ttmp->len=back_tmp->len;
                	memcpy(back_ttmp->val,back_tmp->val,back_tmp->len*sizeof(bit64_t));
                    back_ttmp->num=back_tmp->num;
                    back_ttmp->IDInfo=(struct ID_type*)malloc((back_tmp->num+1)*sizeof(struct ID_type));
                    for(i=0;i<back_tmp->num;i++){
                        back_ttmp->IDInfo[i].seq_id=back_tmp->IDInfo[i].seq_id;
                        back_ttmp->IDInfo[i].pos=back_tmp->IDInfo[i].pos;
                        back_tmp->IDInfo[i].read_num=back_tmp->IDInfo[i].read_num;
                    }
					HASH_ADD(hh,BackInfo_New,val,keylen2,back_ttmp);
					back_ttmp=NULL;
                }
                HASH_FIND(hh,SeedInfo_New,seed_tmp->val,keylen1,seed_temp);
                if(!seed_temp){
					seed_temp=(seed_type*)malloc(sizeof(seed_type)+seed_tmp->len*sizeof(bit64_t));
                    memset(seed_temp,0, sizeof(seed_type)+seed_tmp->len*sizeof(bit64_t));
                    seed_temp->len=seed_tmp->len;
                    memcpy(seed_temp->val,seed_tmp->val,seed_tmp->len*sizeof(bit64_t));
                    seed_temp->num=seed_tmp->num;
                    HASH_ADD(hh,SeedInfo_New,val,keylen1,seed_temp);
					seed_temp=NULL;
                }
                flag++;
                ID_arr=(int*)calloc(back_tmp->num,sizeof(int));index=0;
                ID_arr[index]=back_tmp->IDInfo[0].seq_id;
                i=1;
                while(i<back_tmp->num){
                    if(back_tmp->IDInfo[i].seq_id==ID_arr[index])i++;
                    else{
                        ID_arr[++index]=back_tmp->IDInfo[i].seq_id;i++;
                    }
                }
                index+=1;
                HASH_DEL(BackInfo,back_tmp);free(back_tmp);
                HASH_DEL(SeedInfo,seed_tmp);free(seed_tmp);
            }else{
                HASH_DEL(BackInfo,back_tmp);free(back_tmp);
                HASH_DEL(SeedInfo,seed_tmp);free(seed_tmp);
                HASH_SIZE=HASH_COUNT(SeedInfo);
                continue;
            }
		}
		for(back_tmp=BackInfo;back_tmp!=NULL;back_tmp=(back_type*)(back_tmp->hh.next)){
			keylen2=back_tmp->len*sizeof(bit64_t);
            HASH_FIND(hh,SeedInfo,back_tmp->val,keylen2,seed_temp);
            if(seed_temp){
                Flag_ID=(int*)calloc(back_tmp->num+1,sizeof(int));
                for(j=0;j<back_tmp->num;j++)Flag_ID[j]=1;
                New_Num=0;
                for(j=0;j<back_tmp->num;j++){
                    Find_P=(int*)bsearch(&back_tmp->IDInfo[j].seq_id,ID_arr,index,sizeof(int),compare3);
                    if(!Find_P)New_Num+=back_tmp->IDInfo[j].read_num;
                    else Flag_ID[j]=0;
                }
                if(New_Num>0){
                    int tmp_num=back_tmp->num;
                    for(i=j=0;i<back_tmp->num;i++){
                        if(Flag_ID[i]==1){
                            back_tmp->IDInfo[j].seq_id=back_tmp->IDInfo[i].seq_id;
                            back_tmp->IDInfo[j].pos=back_tmp->IDInfo[i].pos;
                            back_tmp->IDInfo[j].read_num=back_tmp->IDInfo[i].read_num;
                            j++;
                        }else tmp_num--;
                    }
                    back_tmp->num=tmp_num;
                    seed_temp->num=New_Num;
                }else{
                    HASH_DEL(SeedInfo,seed_temp);free(seed_temp);
                    HASH_DEL(BackInfo,back_tmp);free(back_tmp);
                }
                if(Flag_ID){free(Flag_ID);Flag_ID=NULL;}
            }
        }
        if(ID_arr){free(ID_arr);ID_arr=NULL;}
        HASH_SIZE=HASH_COUNT(SeedInfo);
    }
	/*Free memory*/
	back_type *tmp;

	HASH_ITER(hh,SeedInfo,seed_tmp,temp) {
      HASH_DEL(SeedInfo,seed_tmp);
      free(seed_tmp);
    }
	HASH_ITER(hh,BackInfo,back_tmp,tmp) {
      HASH_DEL(BackInfo,back_tmp);
      free(back_tmp);
    }
	/*Output results*/
    printf("Start outputing file\n");
    Output_Seed(Readnum,Uniq_Seed_Num);

	/*Free memory*/
	HASH_ITER(hh,SeedInfo_New,seed_tmp,temp) {
      HASH_DEL(SeedInfo_New,seed_tmp);
      free(seed_tmp);
    }
	HASH_ITER(hh,BackInfo_New,back_tmp,tmp) {
      HASH_DEL(BackInfo_New,back_tmp);
      free(back_tmp);
    }
}
/*Deal input file with fastq.gz*/
void Deal_GZfile()
{
	gzFile gzfp;
	char line[1000];
    int line_num,index,t_id,t_num;
	
	if((gzfp=gzopen(InFileName,"r"))==NULL){
        fprintf(stderr,"Can't open the input file!\n");
        exit(EXIT_FAILURE);
    }
	line_num=0;
    while(gzgets(gzfp,line,sizeof(line)))line_num++;
	SeqInfo=(struct seq_t*)calloc((int)(line_num/2)+1,sizeof(struct seq_t));
	gzseek(gzfp,0,SEEK_SET);
	index=0;
	while(gzgets(gzfp,line,sizeof(line))){
		chomp(line);
		if(line[0]=='>'){
			t_id=Get_Seq_ID(line);t_num=Get_Seq_Num(line);
		}else{
			if(memchr(line,'N',strlen(line)))continue;
			SeqInfo[index].id=t_id;
			SeqInfo[index].num=t_num;
			SeqInfo[index].seq=(char*)calloc(strlen(line)+1,sizeof(char));
			strcpy(SeqInfo[index].seq,line);index++;
		}
	}
	gzclose(gzfp);
	/*Fetch Seed from each read*/
	qsort(SeqInfo,index,sizeof(struct seq_t),compare1);
	Cluster_Seed(index);

	/*Free Memory*/
	int i;
	if(SeqInfo){
		for(i=0;i<index;i++){
			free(SeqInfo[i].seq);SeqInfo[i].seq=NULL;
		}
		free(SeqInfo);SeqInfo=NULL;
	}
}
/*Deal input file with fastq format*/
void Deal_File()
{
	FILE *fp;
	char line[1000];
    int line_num,index,t_id,t_num;

	if((fp=fopen(InFileName,"r"))==NULL){
        fprintf(stderr,"Can't open input file!\n");
        exit(EXIT_FAILURE);
    }
	line_num=0;
    while(fgets(line,sizeof(line),fp))line_num++;
    SeqInfo=(struct seq_t*)calloc((int)(line_num/2)+1,sizeof(struct seq_t));
    fseek(fp,0,SEEK_SET);
	index=0;
    while(fgets(line,sizeof(line),fp)){
		chomp(line);
		if(line[0]=='>'){
			t_id=Get_Seq_ID(line);t_num=Get_Seq_Num(line);
        }else{
			if(memchr(line,'N',strlen(line)))continue;
			SeqInfo[index].id=t_id;
            SeqInfo[index].num=t_num;
            SeqInfo[index].seq=(char*)calloc(strlen(line)+1,sizeof(char));
            strcpy(SeqInfo[index].seq,line);
			index++;
        }
    }
    fclose(fp);
	/*Fetch Seed from each read*/
	qsort(SeqInfo,index,sizeof(struct seq_t),compare1);
    Cluster_Seed(index);

	/*Free Memory*/
    int i;
    if(SeqInfo){
        for(i=0;i<index;i++){
            free(SeqInfo[i].seq);SeqInfo[i].seq=NULL;
        }
        free(SeqInfo);SeqInfo=NULL;
    }
}
int main(int argc,char **argv)
{
    if(argc<7)Usage();

	/*Initialize parameters*/
	SeedLen=40;MaxSeedNum=400;MinSeedNum=2;
	FMaxSeedNum=10*MaxSeedNum;MinSeqNum=2;VMaskNum=10;JMaskNum=5;
    optInit(argc,argv);

	/*Calculate the number of bit64_t seed corresponding to char seed*/
	if(SeedLen%MAX_LEN==0)ISeedLen=(int)SeedLen/MAX_LEN;
	else ISeedLen=(int)SeedLen/MAX_LEN+1;

	/*Open Read Files*/
    char *file_flag;

    file_flag=(char*)calloc(3,sizeof(char));
    strncpy(file_flag,InFileName+strlen(InFileName)-2,2);
    if(!strcmp(file_flag,"gz"))Deal_GZfile();
    else Deal_File();
    if(file_flag){free(file_flag);file_flag=NULL;}

	exit(EXIT_SUCCESS);
}
