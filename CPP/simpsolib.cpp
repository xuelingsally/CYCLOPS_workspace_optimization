/*
 * FILE: simpsolib.cpp, v.1.7.1, 4/28/2014
 * Author: Tomas V. Arredondo
 *
 * SimPSOLib: A simple yet flexible PSO implementation in C++.
 *
 * DISCLAIMER: No liability is assumed by the author for any use made
 * of this program.
 * DISTRIBUTION: Any use may be made of this program, as long as the
 * clear acknowledgment is made to the author in code and runtime executables
 */
#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <stdarg.h>
#include <stdlib.h>
#include "simpsolib.h"

using namespace std;
using namespace simpsolib;
using namespace simtstlib;
using namespace Eigen;

// create writeFile variable
ofstream writeFile("out.txt");

int function_cnt = 0;  // function eval counter

// seed random number generator for ran2
long semran2 = -1000;

bool simpsolib::Population::evaluate()
{
    double fn_value;
    bool success = false;

    // evaluate fitness
    for (std::vector<Organism*>::iterator it_pool = pool.begin(); it_pool != pool.end(); ++it_pool)
    {
        // This is the function for which the population is trying to find the max value.
        //fn_value = evaluator.evaluate((*it_pool)->position);

        fn_value = evaluator.evaluate((*it_pool)->position, evaluator.Input);
        function_cnt++;
        (*it_pool)->value=fn_value;

        //std::cout << fn_value << std::endl;

        if (fn_value > (*it_pool)->best_value)
        {
            (*it_pool)->best_position=(*it_pool)->position;
            (*it_pool)->best_value=fn_value;
        }

        if (fn_value > overall_best_value)
        {
            overall_best_position=(*it_pool)->position;

            overall_best_value=fn_value;

            //pop_leader = it_pool; //Record the population leader
            pop_leader_index = it_pool - pool.begin();


            success = true;
            
            /*
            cout << "Leading Position is: " << endl;
            for (int i=0; i<num_dims; i++)
                cout <<  (*pool[pop_leader_index]).position[i] << ", ";

            cout << endl;
            cout << "Leading Value is: " << (*pool[pop_leader_index]).value << endl; 
            cout << "Leading Index is: " << pop_leader_index << endl;
            */
        }
    }
    return success;
}


void simpsolib::Population_data::evaluate_population_info(Population *pop)
{

    int pop_size=pop->getSize();
    clear_pop_data();
    int cnt=0;

    min_value=(*(pop->pool[0])).value;
    min_index=0;

    for (std::vector<Organism*>::iterator it_pool = pop->pool.begin(); it_pool != pop->pool.end(); ++it_pool)

    {
        sum_values+=(*it_pool)->value;

        if ((*it_pool)->value > max_value)
        {
            max_value=(*it_pool)->value;
            // max_organism=(*(pop->pool[i]));
            max_index=cnt;
        }
        if ((*it_pool)->value < min_value)
        {
            min_value=(*it_pool)->value;
            // min_organism=(*(pop->pool[i]));
            min_index=cnt;
        }
        cnt++;
    }
    if (sum_values)
        avg_value=sum_values/pop_size;
    else
        avg_value=0;

}

void simpsolib::Population_data::display_population_stats()
{
    vector<int>::const_iterator it;
#if 0
    for (it = max_organism.genes.begin(); it != max_organism.genes.end(); it++)
    {
        std::cout << "fittest member ["<< max_index << "]";
        std::cout << *it << " ";
        std::cout << std::endl;
    }
#endif

    std::cout << "fittest member ["<< max_index << "]" << endl;
    std::cout << "fittest member value: "<< max_value << endl;

#if 0
    for (it = min_organism.genes.begin(); it != min_organism.genes.end(); it++)
    {
        std::cout << "least fit member ["<< min_index << "]:";
        std::cout << *it << " ";
        std::cout << std::endl;
    }

    std::cout << "least fit member value: "<< min_value;
    std::cout << std::endl;
#endif

    std::cout << "pop avg value: "<< avg_value;
    std::cout << std::endl;
}

void simpsolib::Population::destroy()
{

    // destroy population
    for ( vector<Organism*>::iterator i = pool.begin(); i != pool.end(); ++i )
    {
        delete *i;
    }
    // empty the container
    pool.clear();
}

void simpsolib::Population::create()
{
    overall_best_value=-2.0;
    pool.resize(population_size);

    for (int i=0; i < population_size; i++)
        pool[i]=new Organism(num_dims);

    // initialize population
    for (int i = 0; i < population_size; i++)
    {
        for (int j = 0; j < num_dims; j++)
        {
            (*pool[i]).position[j]=ran2((double)evaluator.lower_range[j],(double)evaluator.upper_range[j]);
            (*pool[i]).velocity[j]=0;
        }
    }
/*
    (*pool[0]).position[0] = 4.0349;
    (*pool[0]).position[1] = 6.0938;
    (*pool[0]).position[2] = 1.8321;
    (*pool[0]).position[3] = 4.3546;
    (*pool[0]).position[4] = 6.2832;
    (*pool[0]).position[5] = 2.0008;

    (*pool[0]).position[6] = -7.2716;
    (*pool[0]).position[7] = -11.3936;
    (*pool[0]).position[8] = -21.7060;
    (*pool[0]).position[9] = 11.5591;
    (*pool[0]).position[10] = 16.7089;
    (*pool[0]).position[11] = 19.7650;

    (*pool[0]).position[12] = -85.1078;
    (*pool[0]).position[13] = -4.9541;
    (*pool[0]).position[14] = 70.0;
*/
    pop_leader_index = 0;
    rand_resample_count = 0;
}

void simpsolib::Population::display()
{
    for (int i = 0; i < population_size; i++)
    {
        std::cout << "member ["<< i << "]";
        // std::cout << pool[i]; TODO
        std::cout << std::endl;
        std::cout << "value: "<< (*pool[i]).value << endl;
        // std::cout << "double_genes: "<< pool[i].double_genes << endl;
    }
}

#if 0
void simpsolib::Population::display_fittest()
{
    for (int i = 0; i < population_size; i++)
    {
        std::cout << "member ["<< i << "]";
        std::cout << pool[i];
        std::cout << std::endl;
        std::cout << "value: "<< pool[i].value << endl;
        // std::cout << "double_genes: "<< pool[i].double_genes << endl;
    }
}
#endif

double simpsolib::EvalFN::evaluate(vector<double> position)
{
    int tmp_num_dims=position.size();
    double tmp_position[tmp_num_dims];

    for (int i=0; i< tmp_num_dims; i++)
        tmp_position[i]=position[i];

    return ((*eval_fn)(num_parms,tmp_position));
}

double simpsolib::EvalFN::evaluate(vector<double> position, cyclops::fnInputs Input)
{
    int tmp_num_dims=position.size();
    Matrix<double, 15, 1> eaB;

    for (int i=0; i< tmp_num_dims; i++)
        eaB(i,0) = position[i];

    double temp_result = cyclops::objective_function(eaB, Input.W, Input.f_ee_vec, Input.phi_min,
                              Input.phi_max, Input.t_min, Input.t_max,
                              Input.taskspace, Input.radius_tool, Input.radius_scaffold,
                              Input.length_scaffold);
    //std::cout << temp_result << std::endl;
    return temp_result;

}

#if 0
std::cout << "---------- Unsorted GA Population (press enter) ------" << std::endl << flush;
//   std::cin.get(temp);
for (int i = 0; i < population_size; i++)
{
    std::cout << "member ["<< i << "]" << endl;
    std::cout << "member index:" << pool[i].index << endl;
    std::cout << pool[i];
    std::cout << std::endl;
    std::cout << "value: "<< pool[i].value << endl;
    //std::cout << "double_genes: "<< pool[i].double_genes << endl;
}
std::cout << flush;
char temp2;
std::cin.get(temp2);
#endif

void simpsolib::Population::update_vel()
{
    double r_p;
    double r_g;

    // evaluate value
    for (std::vector<Organism*>::iterator it_pool = pool.begin(); it_pool != pool.end(); ++it_pool)
    {
        r_p=ran2(&(semran2));
        r_g=ran2(&(semran2));

        for (int i=0; i< num_dims; i++) // Shi, Eberhart (1998, 2001)
        {
            //std::cout << overall_best_position[i] << std::endl;
            (*it_pool)->velocity[i] = omega*((*it_pool)->velocity[i]) + phi_p*r_p*((*it_pool)->best_position[i] - (*it_pool)->position[i]) + phi_g*r_g*(overall_best_position[i] - (*it_pool)->position[i]);
            
            // Limit velocity in each dimension to full dynamic range in search space (Simplifying PSO, Pedersen 2009)
            /*if (fabs((*it_pool)->velocity[i]) > (evaluator.upper_range[i] - evaluator.lower_range[i]))
            {
                (*it_pool)->velocity[i]=ran2(&(semran2))*(evaluator.upper_range[i] - evaluator.lower_range[i]);
            }*/
            if (fabs((*it_pool)->velocity[i]) > evaluator.max_velocity_vector[i])
            {
                (*it_pool)->velocity[i]=evaluator.max_velocity_vector[i];
            }
        }
    }

}

void simpsolib::Population::update_pos()
{
    // evaluate value
    for (std::vector<Organism*>::iterator it_pool = pool.begin(); it_pool != pool.end(); ++it_pool)
    {
        for (int i=0; i< num_dims; i++) // Shi, Eberhart (1998, 2001)
        {
            (*it_pool)->position[i] = (*it_pool)->position[i] + (*it_pool)->velocity[i] ;
   
            // Limit velocity in each dimension to full dynamic range in search space (Simplifying PSO, Pedersen 2009)
            if ((*it_pool)->position[i] > (evaluator.upper_range[i]))
            {
                (*it_pool)->position[i]=evaluator.upper_range[i];
            }
            else if ((*it_pool)->position[i] < (evaluator.lower_range[i]))
            {
                (*it_pool)->position[i]=evaluator.lower_range[i];
            }

            //if (it_pool == pool.begin())
            //    std::cout << (*it_pool)->position[i] << ", ";
        }
    }
    //std::cout << std::endl;
}

void simpsolib::Population::rand_resample()
{
    for (std::vector<Organism*>::iterator it_pool = pool.begin(); it_pool != pool.end(); ++it_pool)
    {
        double random_num = ran2(&(semran2));
        double likelihood = 1.5 + (*it_pool)->value;
        if (random_num >= likelihood)
        {
            rand_resample_count++;
            for(int i=0; i< num_dims; i++)
            {
                (*it_pool)->position[i] = ran2((double)evaluator.lower_range[i],(double)evaluator.upper_range[i]);
                (*it_pool)->velocity[i] = 0;
            }
        }
    }
}
//-----------------------------------------------------------------------------
// Pattern Search
//-----------------------------------------------------------------------------
void simpsolib::Population::initpatternsearch()
{
    // Determine the mesh size as the minimum of the search space x the search factor
    initial_search_factor = 0.2;
    double temp_min = 999999999;
    for (int i=0; i<num_dims; i++)
    {
        double diff = evaluator.upper_range[i] - evaluator.lower_range[i];
        if (diff < temp_min)
        {
            temp_min = diff;
        }
    }

    if (temp_min <= 1e-6 * 2.0)
    {
        temp_min = 1e-6 * 2.0;
    }

    mesh_size = temp_min * initial_search_factor;

    //mesh_size = 1;

    // create the search vector for pattern search
    
    SearchDirVec = vector<double> (num_dims); 
    
    for (int i=0; i<num_dims; i++)
    {
        //SearchDirVec[i] = (evaluator.upper_range[i] - evaluator.lower_range[i]) * 0.05;
        SearchDirVec[i] = 1.0;
    }
}

void simpsolib::Population::patternsearch()
{

    //cout <<  "Pattern Search! " << endl;
    double fn_value;
    bool success = false;
    // the Mesh set to evaluate is the positive and negative of each direction in the search space
    /*
    cout << "Leading Index is: " << pop_leader_index << endl;
    cout << "Current Leading Positon is: " << endl;
    for (int j=0; j<num_dims; j++)
    {
        cout << (*pool[pop_leader_index]).best_position[j] << ", ";
    }
    cout << endl << endl; */

    for(int i=0; i<num_dims; i++)
    {
        // Create the +ve search vector, d_plus
        vector<double> d_plus = (*pool[pop_leader_index]).best_position;
        d_plus[i] = d_plus[i] + mesh_size * SearchDirVec[i];

        if (d_plus[i] > evaluator.upper_range[i])
        {
            d_plus[i] = evaluator.upper_range[i];
        }
        else if (d_plus[i] < evaluator.lower_range[i])
        {
            d_plus[i] = evaluator.lower_range[i];
        }


        fn_value = evaluator.evaluate(d_plus, evaluator.Input);
        function_cnt++;
/*
        for (int j=0; j<num_dims; j++)
        {
            cout << d_plus[j] << ", ";
        }
        cout << endl << endl; */


        if (fn_value > (*pool[pop_leader_index]).best_value)
        {
            (*pool[pop_leader_index]).position = d_plus;
            (*pool[pop_leader_index]).best_position = d_plus;
            (*pool[pop_leader_index]).best_value = fn_value;
            (*pool[pop_leader_index]).value = fn_value;

            //since this is the pop leader, it is also the best value
            overall_best_position=d_plus;
            overall_best_value=fn_value;

            success = true;
            break;
        }



        // Create the -ve search vector, d_minus
        vector<double> d_minus = (*pool[pop_leader_index]).best_position;
        d_minus[i] = d_minus[i] - mesh_size * SearchDirVec[i];

        if (d_minus[i] > evaluator.upper_range[i])
        {
            d_minus[i] = evaluator.upper_range[i];
        }
        else if (d_minus[i] < evaluator.lower_range[i])
        {
            d_minus[i] = evaluator.lower_range[i];
        }


        fn_value = evaluator.evaluate(d_minus, evaluator.Input);
        function_cnt++;

        if (fn_value > (*pool[pop_leader_index]).best_value)
        {
            (*pool[pop_leader_index]).position = d_minus;
            (*pool[pop_leader_index]).best_position = d_minus;
            (*pool[pop_leader_index]).best_value = fn_value;
            (*pool[pop_leader_index]).value = fn_value;

            //since this is the pop leader, it is also the best value
            overall_best_position = d_minus;
            overall_best_value = fn_value;

            success = true;
            break;
        }
    }


    if (!success)
    {
        // reduce the mesh size if the poll step was not successful
        if (mesh_size > 10e-6)
        {
            mesh_size = mesh_size * 0.5;
        }
        //cout << "Patternsearch Failed, MeshSize: " << mesh_size << endl;
        //writeFile << "  Patternsearch Failed, MeshSize: " << mesh_size << endl;
    }
    else
    {
        // increase the mesh size
        //if (mesh_size < 4)
        //{
        mesh_size = mesh_size * 2.0;
        //}
        /*
        cout << "New Leading Positon is: " << endl;
        for (int j=0; j<num_dims; j++)
        {
            cout << (*pool[pop_leader_index]).position[j] << ", ";
        }
        cout << endl << endl;
        */

        //cout << "Patternsearch Successful, MeshSize: " << mesh_size << endl;
        //writeFile << "  Patternsearch Successful, MeshSize: " << mesh_size << endl;
    }

}
//-----------------------------------------------------------------------------
// Optimization Methods
//-----------------------------------------------------------------------------
int simpsolib::run_pso(EvalFN eval, int number_runs, int pso_pop_size, int pso_number_iters,
                       float phi_p, float phi_g, double omega_initial, double omega_final, bool rand_particle_upd_flag)
{
    clock_t start,end,diff=0;
    start=clock();

    // Run params
    int nRun=0;
    double nRunsAvgFitness=0;
    double nRunsMaxFitness=-100000;
    double nRunsMinFitness=999999999;
    double nRunFitness[number_runs];
    double nRunsAvgFitnessVariance=0;
    double nRunsAvgFitnessSTD=0;
    function_cnt=0;

    for (int i=0; i < number_runs ; i++)
        nRunFitness[i]=0;

    // Set the population parameters: number of dimensions and other parms
    Population pop(eval.num_parms);

    // Set pop parameters
    pop.setEvalFN(eval);
    Population_data pop_info;
    pop.setSize(pso_pop_size);
    pop.setNumIters(pso_number_iters);

    srand(clock());


    for (nRun=0; nRun < number_runs; nRun++)
    {

        //    std::cout << "---------- RUN NUMBER" << nRun << std::endl;
        //    std::cout << "---------- BEGIN OPTIMIZATION ----------" << std::endl;

        pop.create(); // instantiate (new) the population
        pop.setPhiP(phi_p);
        pop.setPhiG(phi_g);
        pop.setOmega(omega_initial);
        pop.omega_initial = omega_initial;
        pop.omega_final = omega_final;
        pop.setRandPartUpdFlag(rand_particle_upd_flag);
        pop.initpatternsearch(); //initialise pattern search parameters

        pop.evaluate();
        pop_info.evaluate_population_info(&pop);


        // std::cout << "---------- Initial Population Stats (press enter) -------" << std::endl << flush;
        // std::cin.get(temp);
        // pop_info.display_population_stats();

        for (int i=1; i < pop.getNumIters() ; i++)
        {
            
            std::cout << "Iteration: " << i << "  Leader: " << pop.pop_leader_index << "  ObjVal: " << pop.getBestVal() << "  meshsize: " << pop.mesh_size << "  RanResampleCount: " << pop.rand_resample_count << std::endl;
            writeFile << "Iteration: " << i << "  Leader: " << pop.pop_leader_index << "  ObjVal: " << pop.getBestVal() << "  meshsize: " << pop.mesh_size << "  RanResampleCount: " << pop.rand_resample_count << std::endl;
            
            double omega_temp = omega_initial - ((omega_initial-omega_final) * i/pop.getNumIters());
            pop.setOmega(omega_temp);
            pop.update_vel();
            pop.update_pos();

            bool temp_success = pop.evaluate();

            // perform pattern search of particle swarm failed.
            if (!temp_success)
            {
                pop.patternsearch();
            }

            pop.rand_resample();
            pop_info.evaluate_population_info(&pop);


            if (i%10==0)
            {
                vector<double> temp_position = (pop.getBestPos());
                writeFile << std::endl << "eaB = [";
                for (int j=0; j < pop.getNumDims(); j++)
                {
                    writeFile << temp_position[j];
                    if (j < pop.getNumDims()-1)
                        writeFile << "; ";
                }
                writeFile << "];" << std::endl << std::endl;
            }
            //pop_info.display_population_stats();

            // std::cout << "iteration: "<< i << "-- Press enter to continue --" << std::endl << flush;
            //std::cin.get(temp);
        }

        // Do a final pattern search till. with mesh size as the limiting factor.

        int ps_counter = 0;
        while (pop.mesh_size > 1e-4)
        {
            ps_counter++;
            pop.patternsearch();
            if (ps_counter%20 == 0)
            {
                cout << "PS iteration: " << ps_counter << "  ObjVal: " << pop.getBestVal() << std::endl;
            }
        }

        //std::cout << "---------- Final Population (press enter) ------" << std::endl << flush;
        //std::cin.get(temp);
        //pop.display();
        //std::cout << flush;

        // std::cout << "---------- Final Population Stats (press enter) -------" << std::endl << flush;
        // pop_info.display_population_stats();
        // std::cin.get(temp);
        // std::cout << "---------- END OPTIMIZATION ----------" << std::endl;

        // Add up numbers for each Run
        if (pop_info.max_value > nRunsMaxFitness)
        {
            nRunsMaxFitness = pop_info.max_value;
        }
        if (pop_info.min_value < nRunsMinFitness)
        {
            nRunsMinFitness = pop_info.min_value;
        }

        nRunFitness[nRun]=pop_info.avg_value;
        nRunsAvgFitness+=pop_info.avg_value;
        
        std::cout << "Max Value: " << pop.getBestVal() << std::endl;

        vector<double> temp_position = (pop.getBestPos());
        std::cout << "eaB = [";
        writeFile << "eaB = [";
        for (unsigned int i=0; i < temp_position.size(); i++)
        {
            std::cout << temp_position[i];
            writeFile << temp_position[i];
            if (i < temp_position.size()-1)
            {
                std::cout << "; ";
                writeFile << "; ";
            }
        }

        std::cout << "];" << std::endl;
        writeFile << "];" << std::endl << std::endl;

        pop.destroy();  // del the population
    }

    end=clock();
    diff=end-start;

    // Calculate RUN results
    nRunsAvgFitness=nRunsAvgFitness/number_runs;

    for (int i=0; i<nRun; i++)
    {
        nRunsAvgFitnessVariance+=pow(nRunsAvgFitness-nRunFitness[i],2);
    }
    nRunsAvgFitnessVariance=nRunsAvgFitnessVariance/number_runs;
    nRunsAvgFitnessSTD=pow(nRunsAvgFitnessVariance, 0.5);

    std::cout << "---------- PSO FINAL RESULTS ----------" << std::endl;
    std::cout << "Number of Runs: "  << nRun << std::endl;
    std::cout << "Function: "        << eval.szName << std::endl;
    std::cout << "Runs max value: "<< nRunsMaxFitness << std::endl;
    std::cout << "Runs min value: "<< nRunsMinFitness << std::endl;
    std::cout << "Runs avg value: "<< nRunsAvgFitness << std::endl;
    std::cout << "Runs std dev value: "<< nRunsAvgFitnessSTD << std::endl;
    std::cout << "Simulation time in clocks: " << diff << std::endl;
    std::cout << "Simulation clocks per sec: " << CLOCKS_PER_SEC << std::endl;
    std::cout << "----------- PARAMETERS USED: -----------" << std::endl;
    std::cout << "Pop size:             "   << pop.getSize() << std::endl;
    std::cout << "Number of Dimensions: "     << pop.getNumDims() << std::endl;
    std::cout << "Number of Iterations: "     << pop.getNumIters() << std::endl;
    std::cout << "Number of Function Evals: "     << function_cnt << std::endl;
    std::cout << "---------- END FINAL RUN RESULTS ------------" << std::endl;

    return 0;
}
