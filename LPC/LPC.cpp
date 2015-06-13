#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>

float x[10000] = {0,};
float s[10000] = {0,};
float S[10000] = {0,};
float w[10000] = {0,};
float R[10000] = {0,};
float K[10000] = {0,};

float A[10000] = {0,};
float AA[10000] = {0,};
float SM[10000] = {0,};
int COUNT = 0;

void getAudio(){
	FILE	*fi;
	float	input;
	short	data;

	fi = fopen("Male.raw", "rb");
	for( ; ; ) {
		if(fread(&data, 2, 1, fi) == NULL)
			break;

		input = (float)data;  // DFT에서는 array로 저장
		x[COUNT] = input;
		COUNT++;
		if(COUNT==320)
			break;
	}
	fclose(fi);
}
void dft(){
	float Xr[10000] = {0,};
	float Xi[10000] = {0,};

	int N = 320;

	FILE* fout;
	fout = fopen("lpc.txt", "w");

	//getAudio();

	for(int k=0;k<N;k++){
		for(int n=0;n<N;n++){
			Xr[k] = Xr[k] + s[n]*cos(2*M_PI*(float)k*(float)n/(float)N);
			Xi[k] = Xi[k] + s[n]*sin(2*M_PI*(float)k*(float)n/(float)N);
		}
		S[k] = sqrt(Xr[k] * Xr[k] + Xi[k]*Xi[k]);
		fprintf(fout,"%f\n",S[k]);
	}
	fclose(fout);
}
void dftOfA(){
	float Xr[10000] = {0,};
	float Xi[10000] = {0,};

	int N = 320;

	FILE* fout;
	fout = fopen("lpc.txt", "w");

	//getAudio();

	for(int k=0;k<N;k++){
		for(int n=0;n<N;n++){
			Xr[k] = Xr[k] + A[n]*cos(2*M_PI*(float)k*(float)n/(float)N);
			Xi[k] = Xi[k] + A[n]*sin(2*M_PI*(float)k*(float)n/(float)N);
		}
		AA[k] = sqrt(Xr[k] * Xr[k] + Xi[k]*Xi[k]);
		//fprintf(fout,"%f\n",AA[k]);
	}
	fclose(fout);
}
void setWindow(){
	for(int i=0;i<320;i++){
		w[i]=0.54-0.46*cos(2*M_PI/319*i);
	}
}
void setS(){
	for(int i=0;i<320;i++){
		s[i] = x[i] * w[i];
	}
}
void Autocorrelation(){
	for(int k=0;k<=10;k++){
		for(int n =0;n<320-k;n++){
			R[k] = R[k] + s[n] * s[n+k]; 
		}
	}
}

int LevinsonRecursion()
{
    unsigned int P = 10;
	double Am1[62];

    if(R[0]==0.0) { 
        for(unsigned int i=1; i<=P; i++) 
        {
            K[i]=0.0; 
            A[i]=0.0;
        }}
    else {
        double km,Em1,Em;
        unsigned int k,s,m;
        for (k=0;k<=P;k++){
            A[0]=0;
            Am1[0]=0; }
        A[0]=1;
        Am1[0]=1;
        km=0;
        Em1=R[0];
        for (m=1;m<=P;m++)                    //m=2:N+1
        {
            double err=0.0f;                    //err = 0;
            for (k=1;k<=m-1;k++)            //for k=2:m-1
                err += Am1[k]*R[m-k];        // err = err + am1(k)*R(m-k+1);
            km = (R[m]-err)/Em1;            //km=(R(m)-err)/Em1;
            K[m-1] = -float(km);
            A[m]=(float)km;                        //am(m)=km;
            for (k=1;k<=m-1;k++)            //for k=2:m-1
                A[k]=float(Am1[k]-km*Am1[m-k]);  // am(k)=am1(k)-km*am1(m-k+1);
            Em=(1-km*km)*Em1;                //Em=(1-km*km)*Em1;
            for(s=0;s<=P;s++)                //for s=1:N+1
                Am1[s] = A[s];                // am1(s) = am(s)
            Em1 = Em;                        //Em1 = Em;
        }
    }
    return 0;
}
void setSpectrum(){
	FILE* fout;
	fout = fopen("lpc.txt", "w");
	for(int i=0;i<320;i++){
		SM[i] = 1 / abs(AA[i]);
		fprintf(fout,"%f\n",SM[i]);
	}
	fclose(fout);
}
int main(){
	getAudio();
	setWindow();
	setS();
	dft();
	Autocorrelation();
	LevinsonRecursion();
	dftOfA();
	setSpectrum();
	return 0;
}
