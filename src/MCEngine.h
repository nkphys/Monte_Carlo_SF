#include "ParametersEngine.h"
#include "Coordinates.h"
#include "MFParams.h"
#include "Hamiltonian.h"
#include "Observables.h"

#ifndef MCENGINE_H
#define MCENGINE_H


class MCEngine{
public:
    MCEngine(Parameters& Parameters__, Coordinates& Coordinates__,
             MFParams& MFParams__, Hamiltonian& Hamiltonian__,
             Observables& Observables__)
        : Parameters_(Parameters__),Coordinates_(Coordinates__),
          MFParams_(MFParams__), Hamiltonian_(Hamiltonian__),
          Observables_(Observables__),
          lx_(Parameters_.lx), ly_(Parameters_.ly), ns_(Parameters_.ns),
          orbs_(Parameters_.orbs)
    {

    }

    void RUN_MC();
    double Prob (double muu, double mu_new);

    Parameters &Parameters_;
    Coordinates &Coordinates_;
    MFParams &MFParams_;
    Hamiltonian &Hamiltonian_;
    Observables &Observables_;
    const int lx_, ly_, ns_, orbs_;

};

/*
 * ***********
 *  Functions in Class MCEngine ------
 *  ***********
*/

void MCEngine::RUN_MC(){

    complex<double> zero(0.0,0.0);
    bool Metropolis_Algo = Parameters_.Metropolis_Algorithm;
    bool Heat_Bath_Algo = Parameters_.Heat_Bath_Algorithm;

    int MC_sweeps_used_for_Avg=Parameters_.Last_n_sweeps_for_measurement;
    int Gap_bw_sweeps = Parameters_.Measurement_after_each_m_sweeps;

    double PrevE,CurrE,P_new,P12,muu;
    double Curr_QuantE;
    int x,y,act;
    double saved_Params[2];

    string File_Out_progress;
    string File_Out_theta_phi;

    double temp_=Parameters_.temp_max;



    double Curr_Cluster_CLE;

    //starting with a random guess

    while(temp_>=Parameters_.temp_min){

        cout << "Temperature = " << temp_<<" is being done"<<endl;
        Parameters_.temp=temp_;
        Parameters_.beta=double(11604.0/ temp_);

        for(int ix=0;ix<lx_;ix++){
            for(int iy=0;iy<ly_;iy++){
                Observables_.SiSjQ_Mean_(ix,iy)=zero;
                Observables_.SiSjQ_square_Mean_(ix,iy)=zero;
            }
        }
        Observables_.AVG_Total_Energy=0.0;
        Observables_.AVG_Total_Energy_sqr=0.0;

        MFParams_.etheta_avg.fill(0.0);
        MFParams_.ephi_avg.fill(0.0);



        char temp_char[50];
        sprintf(temp_char,"%.1f",temp_);

        File_Out_progress = "output_Temp" + string(temp_char) + ".txt";
        ofstream file_out_progress(File_Out_progress.c_str());

        File_Out_theta_phi = "ThetaPhi_Temp" + string(temp_char) + ".txt";
        ofstream File_Out_Theta_Phi(File_Out_theta_phi.c_str());



        file_out_progress<< "Total "<<Parameters_.IterMax<<" sweeps are performed."<<endl;
        file_out_progress<<"First "<<Parameters_.IterMax - (Gap_bw_sweeps*(MC_sweeps_used_for_Avg - 1) + MC_sweeps_used_for_Avg)<<
                           " sweeps are used for thermalization and every "<<Gap_bw_sweeps+1<<" in last "<<
                           Gap_bw_sweeps*(MC_sweeps_used_for_Avg - 1) + MC_sweeps_used_for_Avg<<
                           " sweeps are used for measurement."<<endl;
        act=1;




        Parameters_.WindowSize = 0.1; //2f + 0.003f*beta0 ;
        Parameters_.Eav=0.0;
        Parameters_.MCNorm=0;
        Parameters_.Dflag='N'; // flag to calculate only Eigenvalue
        //std::string name="Output/Conf_" + to_string(ltemp) + ".dat";
        //Parameters_.beta = double(11604.0/ (Parameters_.temp +20.0) );
        //cout << "TEMP  " << Parameters_.temp << endl;


        file_out_progress<<"I_MC"<<setw(15)<<"S(0,1)"<<setw(15)<<"S(1,0)"
                        <<setw(15)<<"S(0,Pi)"<<setw(15)<<"S(Pi,0)"<<setw(17)<<"< N_total >"
                       <<setw(15)<<"E_CL"<<setw(15)<<"E_QM"<<setw(15)<<"mu"<< endl;




        /*
        Hamiltonian_.InteractionsCreate();
        Hamiltonian_.Diagonalize(Parameters_.Dflag);
        PrevE=Hamiltonian_.GetCLEnergy();
        Hamiltonian_.copy_eigs(1);
        Parameters_.mus=0.25;
        Parameters_.mus=Hamiltonian_.chemicalpotential(0.25,Parameters_.Fill);
        //P_old=Prob(PrevE,Parameters_.mus);
        Observables_.SiSjFULL();

        */





        int Confs_used=0;
        int measure_start=0;

        for(int count=0;count<Parameters_.IterMax;count++){
            //if (count == 1){
            // Parameters_.beta = double(11604.0/ Parameters_.temp);
            // PrevE = Hamiltonian_.GetCLEnergy();
            //  Hamiltonian_.InteractionsCreate();
            //  Hamiltonian_.Diagonalize(Parameters_.Dflag);
            //  Hamiltonian_.copy_eigs(1);
            //muu = Hamiltonian_.chemicalpotential(Parameters_.mus,Parameters_.Fill);
            //Parameters_.mus = Parameters_.mus*0.4f + muu*0.6f;
            // }

            for(int i=0;i<ns_;i++) {  // For each site

                x=Coordinates_.indx(i);
                y=Coordinates_.indy(i);
                saved_Params[0]=MFParams_.etheta(x,y);
                saved_Params[1]=MFParams_.ephi(x,y);

                MFParams_.FieldThrow(i);
                CurrE = Hamiltonian_.GetCLEnergy();
                Hamiltonian_.InteractionsCreate();
                Hamiltonian_.Diagonalize(Parameters_.Dflag);
                muu=Hamiltonian_.chemicalpotential(Parameters_.mus,Parameters_.Fill);
                Curr_QuantE = Hamiltonian_.E_QM();

                //Ratio of Quantum partition functions
                // i.e <exp(-beta(Hquant_new))>/<exp(-beta(Hquant_old))>
                P_new = Prob(Parameters_.mus, muu);


                //P12 = [ <exp(-beta(Hquant_new))>/<exp(-beta(Hquant_old))> ]*
                //      [exp(-beta*E_classical(New)) / exp(-beta*E_classical(old))]
                //      * [sin(Theta_i(New)) / sin(Theta_i(Old)) ]
                P12 = P_new*exp(-Parameters_.beta*(CurrE-PrevE));
                P12*= (sin(MFParams_.etheta(x,y))/sin(saved_Params[0]));

                //Heat bath algorithm [See page-129 of Prof. Elbio's Book]
                //Heat bath algorithm works for small changes i.e. when P12~1.0
                if (Heat_Bath_Algo){
                    P12 =P12/(1.0+P12);
                }



                //Metropolis Algotithm
                if (Metropolis_Algo){
                    P12=min(1.0,P12);
                }




                /*
       * VON NEUMANN's REJECTING METHOD:
       * Random number < P12 -----> ACCEPT
       * Random number > P12 -----> REJECT
       */


                //ACCEPTED
                if ( MFParams_.random() < P12 ) {
                    Parameters_.AccCount[0]++;
                    PrevE=CurrE;
                    Hamiltonian_.copy_eigs(1);
                    Parameters_.mus = muu;
                    act=1;
                }

                //REJECTED
                else{
                    Parameters_.AccCount[1]++;
                    act=0;
                    MFParams_.etheta(x,y) = saved_Params[0];
                    MFParams_.ephi(x,y)   = saved_Params[1];

                }

                // if ((act == 1) && (count<1000)) {

                //muu = Hamiltonian_.chemicalpotential(Parameters_.mus,Parameters_.Fill);
                //Parameters_.mus = Parameters_.mus*0.999 + muu*0.001;
                //Parameters_.mus = muu;
                //}

            }// site loop

            //      if (act == 1) {
            //       muu = Hamiltonian_.chemicalpotential(Parameters_.mus,Parameters_.Fill);
            //       Parameters_.mus = Parameters_.mus*0.99f + muu*0.01f;
            //      }

            if ( (count%10==0) ) {
                MFParams_.Adjust_MCWindow();
            }

            if(count < (Parameters_.IterMax - (Gap_bw_sweeps*(MC_sweeps_used_for_Avg - 1) + MC_sweeps_used_for_Avg)) ){
                if ( (count%10==0) ) {
                    Observables_.SiSjFULL();
                    file_out_progress << int(1.0*count) <<setw(20)<< Observables_.SiSj(0,1) <<setw(16)<< Observables_.SiSj(1,0)
                                      <<setw(16)<< Observables_.SiSjQ(0,int(lx_/2)).real() <<setw(16)<< Observables_.SiSjQ(int(lx_/2),0).real()
                                     <<setw(16)<< Hamiltonian_.TotalDensity() <<setw(16)<< CurrE
                                    <<setw(16)<< Curr_QuantE<<setw(15)<<muu<< endl;
                }
            }
            //Average and Std. deviation is calculated is done
            else{

                if(measure_start==0){
                    measure_start++;
                    file_out_progress<<"----------Measurement is started----------"<<endl;
                    file_out_progress<<"Avg{S(pi,0)}  Avg{S(0,pi)}  std.dev{S(pi,0)}  std.dev{S(0,pi)}  Avg{E}  std.dev{E}"<<endl;
                }
                int temp_count=count -
                        (Parameters_.IterMax - (Gap_bw_sweeps*(MC_sweeps_used_for_Avg - 1) + MC_sweeps_used_for_Avg));
                int zero_or_not = temp_count % (Gap_bw_sweeps + 1);
                if( zero_or_not==0 ){
                    Confs_used=Confs_used+1;
                    Observables_.SiSjFULL();
                    Observables_.SiSjQ_Average();
                    Observables_.Total_Energy_Average(Curr_QuantE, CurrE);

                    MFParams_.Calculate_Fields_Avg();

                    //double MC_steps_Avg_insitu = (1.0 + 1.0*(count - (Parameters_.IterMax - MC_steps_used_for_Avg)));

                    file_out_progress << int(1.0*count) <<setw(20)<< Observables_.SiSjQ_Mean(int(lx_/2),0).real()/(Confs_used*1.0)
                                      <<setw(16)<<Observables_.SiSjQ_Mean(0,int(lx_/2)).real()/(Confs_used*1.0)
                                     <<setw(16)<<
                                       sqrt(
                                       (( Observables_.SiSjQ_square_Mean(int(lx_/2),0)/(Confs_used*1.0) ) -
                                        ((Observables_.SiSjQ_Mean(int(lx_/2),0)*Observables_.SiSjQ_Mean(int(lx_/2),0) )/(Confs_used*Confs_used*1.0) ) ).real()
                                           )

                                    <<setw(16)<<
                                      sqrt(
                                      (( Observables_.SiSjQ_square_Mean(0,int(lx_/2))/(Confs_used*1.0) ) -
                                       ((Observables_.SiSjQ_Mean(0,int(lx_/2))*Observables_.SiSjQ_Mean(0,int(lx_/2)) )/(Confs_used*Confs_used*1.0) ) ).real()
                                          )
                                      <<setw(16)<<
                                       Observables_.AVG_Total_Energy/(Confs_used*1.0)
                                        <<setw(16)<<
                                       sqrt(  (Observables_.AVG_Total_Energy_sqr/(Confs_used*1.0)) -
                                          ((Observables_.AVG_Total_Energy*Observables_.AVG_Total_Energy)/(Confs_used*Confs_used*1.0))  )
                                   <<endl;

                }

            }


        }// Iter Loop
        file_out_progress << "Total "<<Confs_used<< " configurations were used were measurement"<<endl;

        temp_ = temp_ - Parameters_.d_Temp;

        File_Out_Theta_Phi<<"#x"<<setw(15)<<"y"<<setw(15)<<"Theta_avg(x,y)"<<setw(15)<<"Phi_avg(x,y)"<<endl;
        for(int ix=0;ix<lx_;ix++){
            for(int iy =0;iy<ly_;iy++){
                File_Out_Theta_Phi<<ix<<setw(15)<<iy<<setw(15)<<MFParams_.etheta_avg(ix,iy)/(Confs_used*1.0)<<setw(15)<<MFParams_.ephi_avg(ix,iy)/(Confs_used*1.0)<<endl;
            }}


        MFParams_.Read_classical_DOFs(File_Out_theta_phi);
    }//Temperature loop

} // ---------



double MCEngine::Prob(double muu, double mu_new){

    double P=double(1.0), X,Y,X2;

    for(int i=0;i<2*orbs_*ns_;i++){
        X = Parameters_.beta*(mu_new - Hamiltonian_.eigs_[i]);
        Y = Parameters_.beta*(muu - Hamiltonian_.eigs_saved_[i]);
        X2 = exp(-X);
        if ( X < 0.0 ) {
            P=P*((1.0+exp(X))/(1.0+exp(Y)));
        }
        else {
            P=P*((1.0+X2)/(X2+exp(Y-X)));
        }
    }

    return P;

} // ---------







#endif // MCENGINE_H
