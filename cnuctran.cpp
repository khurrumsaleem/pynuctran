// cnuctran.cpp : This file contains all necessary routines and functions that enables the
//                calculation of final nuclide concentrations after depletion processes. 
//                Program execution begins and ends there.
//                This code is prepared by M. R. Omar, Universiti Sains Malaysia.
//

#include <iostream>
#include <iomanip>
#include <mpfr.h>
#include <mpir.h>
#include <mpreal.h>
#include <map>
#include <unordered_map>
#include <vector>
#include <fstream>
#include "pugixml.hpp"
#include <thread>
#include <chrono>
#include <filesystem>

using namespace mpfr;
using namespace std;
using namespace pugi;

namespace cnuctran {

    /*
        REUSABLE HIGH PRECISION CONSTANTS.
        __two__ = 2
        __one__ = 1
        __neg__ = -1
        __zer__ = 0
        __eps__ is the epsilon value. Any value falls below __eps__ is assumed to be zero.

        REUSABLE INTEGER CONSTANTS.
        __dps__ is the mpreal arithmetic precision.
        __dop__ is the decimal places of any printed mpreal numbers.
        __nop__ is an integer specifying no product.
        __vbs__ is the vervosity level; 0 (none), 1 (minimal), 2 (comprehensive)
    */
    const mpreal __two__ = mpreal("2.0", digits2bits(500));
    const mpreal __one__ = mpreal("1.0", digits2bits(500));
    const mpreal __neg__ = mpreal("-1.0", digits2bits(500));
    const mpreal __zer__ = mpreal("0.0", digits2bits(500));
    const mpreal __eps__ = mpreal("1e-200", digits2bits(500));
    const int    __dps__ = 200;
    const int    __dop__ = 16;
    const int    __npr__ = 8;
    const int    __nop__ = -1;
    int          __vbs__ = 0;

    /*
        POINTERS FOR COMPRESSED SPARSE ROW MATRIX MULTIPLICATIONS.
    */





    /*
            This class was initially developed to accomodate fast, high-precision sparse
            matrix multiplications and powers. WARNING! This class does not covers all
            matrix operations, it only cover the basic operations used by CNUCTRAN, i.e.
            Multiplication and Powers.

            This class uses the basic C++ map to store data. Sparse matrix
            elements are accessed at an incredible speed via the use of hash table.

            Of course, there is still no known library that provides high-precision sparse
            matrix operations. Therefore, I must endure writing a new specialized class
            handling sparse matrix power to preserve the accuracy.

            SPARSE STORAGE. Only the non-zero elements are stored in smatrix.data dictio-
            nary. The keys of smatrix.data are the tuple specifying the position of the
            non-zero elements in the dense matrix version. smatrix.common_column (cc) and
            smatrix.common_rows (cr) are dictionaries that stores the collection (also a dict.)
            of position tuples with common column or row indices, respectively. The keys are
            the common column/row indices.

            SPARSE MULTIPLICATION. Consider sparse matrices A and B. We want to evaluate A*B.
            Here we implement the row-wise multiplication algorithm. Each row can be vectorized
            into multiple concurrent threads.

            For a more comprehensive understanding, consider reading the code below. Good luck!

            SPARSE POWER. Suppose we want to evaluate the power of a sparse matrix, i.e. A^n.
            Let n be a large integer number. A naive method is given by,

            A^n = A x A x A x .... (n times)

            Fortunately, this process can be accelerated using the binary decomposition method,
            for instance,

            let C = A x A (power raised to 2)
            C = C x C     (power raised to 4)
            C = C x C     (power raised to 8)
            :
            :
            until...
            C = C x C     (power raised to n)

            This algorithm has a complexity of O(log n).

    */
    unordered_map<int, mpreal>* rdr;
    unordered_map<int, mpreal>* sdd;
    unordered_map<int, mpreal>* odd;

    void row_operation(int& irow,
        unordered_map<int, unordered_map<int, mpreal>>& sd,
        unordered_map<int, unordered_map<int, mpreal>>& od,
        unordered_map<int, unordered_map<int, mpreal>>& rd)
    {


        rdr = &rd[irow];
        sdd = &sd[irow];

        for (unordered_map<int, mpreal>::iterator it1 = sdd->begin(); it1 != sdd->end(); it1++) {
            int icol = it1->first;
            odd = &od[icol];
            mpreal x = (*sdd)[icol];
            for (unordered_map<int, mpreal>::iterator it2 = odd->begin(); it2 != odd->end(); it2++) {
                int ocol = it2->first;
                (*rdr)[ocol] += x * (*odd)[ocol];
            }
        }
    }

    /*void row_operation_chunk(int istart, int istop,
        unordered_map<int, unordered_map<int, mpreal>>& sd,
        unordered_map<int, unordered_map<int, mpreal>>& od,
        unordered_map<int, unordered_map<int, mpreal>>& rd)
    {

        for (int irow = istart; irow <= istop; irow++)
        {
            rdr = &rd[irow];
            sdd = &sd[irow];

            for (unordered_map<int, mpreal>::iterator it1 = sdd->begin(); it1 != sdd->end(); it1++) {
                int icol = it1->first;
                odd = &od[icol];
                mpreal x = (*sdd)[icol];
                for (unordered_map<int, mpreal>::iterator it2 = odd->begin(); it2 != odd->end(); it2++) {
                    int ocol = it2->first;

                    (*rdr)[ocol] += x * (*odd)[ocol];

                }
            }
        }

    }*/

    class smatrix {

    public:

        pair<int, int> shape;
        unordered_map<int, unordered_map<int, mpreal>> data;

        smatrix(pair<int, int> shape)
        {
            this->shape = shape;
            return;
        }

        smatrix(vector<vector<mpreal>>& A)
        {
            int nrows = A.size();
            int ncols = A[0].size();
            this->shape = pair<int, int>(nrows, ncols);

            for (int i = 0; i < nrows; i++)
            {
                for (int j = 0; j < ncols; j++)
                {
                    mpreal a = A[i][j];
                    if (a > __zer__)
                        this->data[i][j] = a;
                }
            }

            return;
        }

        smatrix copy()
        {
            smatrix r = smatrix(this->shape);
            r.data = this->data;
            return r;
        }

        string to_string(int digits = 6)
        {
            ostringstream s;
            ostringstream tmp;
            s.precision(digits);
            for (int i = 0; i < this->shape.first; i++)
            {
                for (int j = 0; j < this->shape.second; j++)
                {
                    mpreal a = this->data[i][j];
                    if (a > __eps__)
                    {
                        tmp.str("");
                        tmp << "(" << i << ", " << j << ")";
                        s << setw(12) << tmp.str() << "\t\t" << a << endl;
                    }
                }
            }

            return s.str();
        }



        smatrix mul(smatrix& other)
        {
            int sx = this->shape.first;
            int sy = other.shape.second;
            smatrix result = smatrix(pair<int, int>(sx, sy));
            unordered_map<int, unordered_map<int, mpreal>>* sd = &this->data;
            unordered_map<int, unordered_map<int, mpreal>>* od = &other.data;
            unordered_map<int, unordered_map<int, mpreal>>* rd = &result.data;

            for (int row = 0; row < sx; row++) {
                row_operation(row, *sd, *od, *rd);

            }

            //vector<thread> thread_pool;
            //int thread_length = sx / __npr__;
            //for (int ithread = 0; ithread < __npr__; ithread++)
            //{
            //    int istart = ithread * thread_length;
            //    int istop = istart + thread_length - 1;

            //    if (ithread != (__npr__ - 1))
            //    {
            //        thread thread_instance(row_operation_chunk, (istart), (istop), ref(*sd), ref(*od), ref(*rd));
            //        thread_instance.join();
            //        
            //    }
            //    else
            //    {
            //        istop = sx - 1;
            //        thread thread_instance(row_operation_chunk, (istart), (istop), ref(*sd), ref(*od), ref(*rd));
            //        thread_instance.join();
            //        
            //    }

            //
            //}



            return result;
        }

        smatrix operator *(smatrix other)
        {
            smatrix r = smatrix(other.shape);
            r = this->mul(other);
            return r;
        }

        smatrix pow(mpz_t n)
        {

            mpz_t zer; mpz_init(zer); mpz_set_ui(zer, 0);
            mpz_t one; mpz_init(one); mpz_set_ui(one, 1);
            mpz_t two; mpz_init(two); mpz_set_ui(two, 2);
            mpz_t modulo; mpz_init(modulo);

            if (mpz_cmp(n, one) <= 0)
            {
                cout << "fatal-error <cnuctran::smatrix::pow()> Matrix exponent is less than one!" << endl;
                exit(1);
            }


            smatrix r = smatrix(this->shape);
            smatrix y = this->copy();
            bool start_flag = true;
            while (mpz_cmp(n, one) > 0)
            {
                mpz_mod(modulo, n, two);
                if (mpz_cmp(modulo, zer) != 0)
                {
                    if (start_flag == true)
                    {
                        start_flag = false;
                        r = y;
                    }
                    else
                        r = r * y;
                }

                mpz_fdiv_q(n, n, two);
                y = y * y;
            }
            if (start_flag)
                r = y;
            else
                r = r * y;
            return r;

        }

    };

    class solver
    {
    public:
        vector<string> species_names;
        int __I__;
        vector<vector<mpreal>> lambdas;
        vector<vector<vector<int>>> G;
        vector<vector<mpreal>> P;
        vector<vector<mpreal>> fission_yields;

        solver(vector<string> species_names)
        {
            this->species_names = species_names;
            this->__I__ = this->species_names.size();
            for (int i = 0; i < this->__I__; i++)
            {
                this->lambdas.push_back(vector<mpreal>());
                vector<int> tmp1 = { __nop__ }; vector<vector<int>> tmp2; tmp2.push_back(tmp1);
                this->G.push_back(tmp2);
                this->P.push_back(vector<mpreal>());
                this->fission_yields.push_back(vector<mpreal>());
            }
            return;
        }

        void add_removal(int species_index,
            mpreal rate,
            vector<int> products,
            vector<mpreal> fission_yields = vector<mpreal>({}))
        {

            // Prints removal info for debugging purpose.
            if (__vbs__ == 2)
            {
                cout << "<define removal> parent = " << this->species_names[species_index];
                cout << " index = " << species_index;
                cout << " rate = " << rate;
                cout << " products = ";

                int i = 0;
                for (int p : products)
                {
                    if (p != __nop__)
                        cout << "[" << p << "] " << this->species_names[p];
                    else
                    {
                        cout << "[" << p << "] not-tracked";
                        continue;
                    }

                    if (fission_yields.size() != 0)
                        cout << " (fission_yield = " << fission_yields[i++] << ")\t";
                    else
                        cout << "\t";
                }
                cout << endl;
            }









            if (rate < __eps__)
                return;

            this->lambdas[species_index].push_back(rate);
            this->G[species_index].push_back(products);

            if (!fission_yields.empty() && products.size() > 1)
            {
                if (fission_yields.size() >= products.size())
                {
                    vector<mpreal> tmp = vector<mpreal>();
                    for (mpreal y : fission_yields)
                        tmp.push_back(y);
                    this->fission_yields[species_index] = tmp;
                }
                else
                {
                    cout << "fatal-error <cnuctran::solver::add_removal(...)> Insufficient fission yields given for species " <<
                        this->species_names[species_index] << " products.";
                    exit(1);
                }
            }
            else if (fission_yields.empty() && products.size() == 1)
            {
                for (int product : products)
                {
                    // This part is for preparing the transmutation matrix, which is not relevant 
                    // for this C++ version, cnuctran. For CRAM calculation, please use the PyNUCTRAN.
                }
            }
            else
            {
                cout << "fatal-error <cnuctran::solver::add_removal(...)> Invalid removal definition for isotope " <<
                    this->species_names[species_index] << endl;
                cout << "Non-fission events MUST only have ONE daughter product." << endl;
                cout << "Whereas fission events MUST have >1 products to track." << endl;
                exit(1);

            }
        }

        smatrix prepare_transfer_matrix(mpreal dt)
        {
            vector<vector<mpreal>> A;
            for (int i = 0; i < this->__I__; i++)
            {
                vector<mpreal> tmp = vector<mpreal>();
                for (int j = 0; j < this->__I__; j++)
                    tmp.push_back(__zer__);
                A.push_back(tmp);
            }

            for (int i = 0; i < this->__I__; i++)
            {
                int n_events = this->G[i].size();

                mpreal norm = __zer__;

                // Compute the probability of removals...
                vector<mpreal> E = vector<mpreal>();
                for (int l = 1; l < n_events; l++)
                    E.push_back(exp(-this->lambdas[i][l - 1] * dt));

                for (int j = 0; j < n_events; j++)
                {
                    this->P[i].push_back(__one__);
                    for (int l = 1; l < n_events; l++)
                    {
                        mpreal kron = mpreal(to_string((int)(l == j)).c_str());
                        this->P[i][j] *= (kron + pow(__neg__, kron) * E[l - 1]);
                    }
                    norm += this->P[i][j];
                }

                if (norm == __zer__)
                    continue;

                for (int j = 0; j < n_events; j++)
                {
                    this->P[i][j] /= norm;
                    int n_daughters = this->G[i][j].size();
                    for (int l = 0; l < n_daughters; l++)
                    {
                        int k = this->G[i][j][l];
                        mpreal a;
                        if (k != __nop__)
                        {
                            if (n_daughters > 1)
                            {
                                a = this->P[i][j] * this->fission_yields[i][l];
                                if (a > __eps__)
                                    A[k][i] += a;
                            }
                            else
                            {
                                a = this->P[i][j];
                                if (a > __eps__)
                                    A[k][i] += a;
                            }

                            if (A[k][i] < __eps__)
                                A[k][i] == __zer__;



                        }
                    }

                    if (j == 0 && this->P[i][j] > __eps__)
                        A[i][i] += this->P[i][j];

                }
            }
            return smatrix(A);
        }

        map<string, mpreal> solve(map<string, mpreal> w0,
            mpreal dt,
            mpreal t)
        {
            vector<vector<mpreal>> w0_matrix;
            for (int i = 0; i < this->__I__; i++)
                w0_matrix.push_back(vector<mpreal>({ __zer__ }));
            for (int i = 0; i < this->__I__; i++)
                if (w0.count(this->species_names[i]) == 1)
                    w0_matrix[i][0] = w0[this->species_names[i]];
            smatrix converted_w0 = smatrix(w0_matrix);
            mpz_t substeps; mpz_set_str(substeps, floor(t / dt).toString().c_str(), 10);
            if (__vbs__) cout << "t = " << t << "\tdt = " << dt << "\tsubsteps = " << substeps << endl;
            auto t1 = chrono::high_resolution_clock::now();
            smatrix A = this->prepare_transfer_matrix(dt);
            
            smatrix An = A.pow(substeps);
            smatrix w = An * converted_w0;
            auto t2 = chrono::high_resolution_clock::now();
            if (__vbs__) cout << "Done computing concentrations. ";
            if (__vbs__) cout << chrono::duration_cast<chrono::milliseconds>(t2 - t1).count() << "ms." << endl;
            map<string, mpreal> out;
            for (int i = 0; i < this->__I__; i++)
                out[this->species_names[i]] = w.data[i][0];

            return out;

        }
    };

    class simulation
    {

    public:

        static enum errex
        {
            MISSING_SUBSTEP_SIZE = 1,
            MISSING_STEP_SIZE = 2,
            MISSING_SPECIES_NAMES = 3,
            MISSING_W0 = 4,
            MISSING_RXN_RATE = 5,
            NUCLIDES_DATA_LOAD_FAILED = 6,
            XML_READING_ERROR = 7,
            UNEXPECTED_ERROR = 8,
            MISSING_W0_SOURCE = 9

        };
        static void build_chains(solver& s, map<string, map<string, mpreal>>& rxn_rates,
            string xml_data_location)
        {
            vector<string> species_names = s.species_names;
            xml_document file;
            xml_parse_result open_success = file.load_file(xml_data_location.c_str());
            if (!open_success)
            {
                cout << "ERROR <cnuctran::depletion_scheme::build_chains(...)> Fail retrieving data from " << xml_data_location << "." << endl;
                return;
            }

            xml_node root = file.child("depletion");

            for (xml_node species : root.children())
            {
                string species_name = species.attribute("name").value();
                vector<string>::iterator it = std::find(species_names.begin(), species_names.end(), species_name);
                if (it == species_names.end())
                    continue;


                mpreal decay_rate;
                if (species.attribute("half_life"))
                    decay_rate = mpfr::log(__two__) / mpreal(species.attribute("half_life").value());
                else
                    decay_rate = __zer__;

                for (xml_node removal : species.children())
                {
                    if (string(removal.name()) == "decay_type")
                    {
                        mpreal decay_rate_adjusted = mpreal(removal.attribute("branching_ratio").value()) * decay_rate;
                        string parent = species_name;
                        string daughter = removal.attribute("target").value();
                        vector<string>::iterator it_parent = std::find(species_names.begin(), species_names.end(), parent);
                        int parent_id = distance(species_names.begin(), it_parent);
                        vector<string>::iterator it_daughter = std::find(species_names.begin(), species_names.end(), daughter);
                        if (it_daughter != species_names.end())
                        {
                            int daughter_id = distance(species_names.begin(), it_daughter);
                            s.add_removal(parent_id, decay_rate_adjusted, vector<int>({ daughter_id }));
                        }
                        else
                        {
                            s.add_removal(parent_id, decay_rate_adjusted, vector<int>({ __nop__ }));
                        }

                    }

                    if (rxn_rates.size() != 0)
                    {
                        if (rxn_rates.count(species_name))
                        {
                            if (string(removal.name()) == "reaction_type" && removal.attribute("target"))
                            {
                                string parent = species_name;
                                vector<string>::iterator it_parent = std::find(species_names.begin(), species_names.end(), parent);
                                int parent_id = distance(species_names.begin(), it_parent);
                                if (rxn_rates[parent].count(removal.attribute("type").value()) &&
                                    removal.attribute("type").value() != "fission")
                                {
                                    string daughter = removal.attribute("target").value();
                                    mpreal removal_rate = mpreal(rxn_rates[parent][removal.attribute("type").value()]);
                                    vector<string>::iterator it_daughter = std::find(species_names.begin(), species_names.end(), daughter);
                                    if (it_daughter != species_names.end())
                                    {
                                        int daughter_id = distance(species_names.begin(), it_daughter);
                                        s.add_removal(parent_id, removal_rate, vector<int>({ daughter_id }));
                                    }
                                    else
                                    {
                                        s.add_removal(parent_id, removal_rate, vector<int>({ __nop__ }));
                                    }
                                }
                            }

                            if (string(removal.name()) == "neutron_fission_yields")
                            {
                                string parent = species_name;
                                vector<string>::iterator it_parent = std::find(species_names.begin(), species_names.end(), parent);
                                int parent_id = distance(species_names.begin(), it_parent);
                                mpreal energy = __zer__;
                                vector<string> products;
                                vector<mpreal> yields;
                                if (rxn_rates[parent].count("fission") != 0)
                                {
                                    for (xml_node data : removal.children())
                                    {
                                        if (string(data.name()) == "energies")
                                        {
                                            stringstream ss(data.child_value()); string token;
                                            vector<string> energies;
                                            while (getline(ss, token, ' '))
                                                if (token != "")
                                                    energies.push_back(token);
                                            energy = energies[0];
                                        }

                                        if (string(data.name()) == "fission_yields")
                                        {
                                            if (mpreal(data.attribute("energy").value()) == energy)
                                            {
                                                for (xml_node param : data.children())
                                                {
                                                    if (string(param.name()) == "products")
                                                    {
                                                        stringstream ss(param.child_value()); string token;
                                                        while (getline(ss, token, ' '))
                                                            if (token != "")
                                                                products.push_back(token);
                                                    }

                                                    if (string(param.name()) == "data")
                                                    {
                                                        stringstream ss(param.child_value()); string token;
                                                        while (getline(ss, token, ' '))
                                                            if (token != "")
                                                                yields.push_back(mpreal(token));
                                                    }

                                                }

                                                mpreal total_fission_rate = rxn_rates[parent]["fission"];
                                                vector<mpreal> yields_to_add;
                                                vector<int> daughters_id_to_add;
                                                vector<string>::iterator it_product;
                                                for (string product : products)
                                                {

                                                    it_product = find(species_names.begin(), species_names.end(), product);
                                                    if (it_product != species_names.end())
                                                    {
                                                        int product_id = distance(species_names.begin(), it_product);
                                                        daughters_id_to_add.push_back(product_id);
                                                        it_product = find(products.begin(), products.end(), product);
                                                        product_id = distance(products.begin(), it_product);
                                                        yields_to_add.push_back(yields[product_id]);
                                                    }
                                                }
                                                it_parent = find(species_names.begin(), species_names.end(), parent);
                                                parent_id = distance(species_names.begin(), it_parent);

                                                s.add_removal(parent_id, total_fission_rate, daughters_id_to_add, yields_to_add);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            return;
        }

        static vector<string> get_nuclide_names(string xml_data_location, int AMin = -1, int AMax = -1)
        {
            xml_document file;
            xml_parse_result open_success = file.load_file(xml_data_location.c_str());
            if (!open_success)
            {
                cout << "ERROR <cnuctran::depletion_scheme::get_nuclide_names(...)> Fail retrieving data from " << xml_data_location << "." << endl;
                return vector<string>();
            }


            xml_node root = file.child("depletion");

            vector<string> species_names;

            for (xml_node species : root.children())
            {
                string name = species.attribute("name").value();
                stringstream ss(name); string token;
                getline(ss, token, '_');
                string x = "";
                for (char c : token)
                    if (isdigit(c)) x += c;
                if (AMin == AMax == -1)
                    species_names.push_back(name);
                else
                {
                    int A = stoi(x);
                    if (A >= AMin && A <= AMax)
                        species_names.push_back(name);
                }
            }
            return species_names;
        }


  /** Runs the simulation according to the given XML input file.
 *
 *  @param xml_input: A string specifying the location of the XML input file.
 *  @return void
 *
 */
        static void from_input(string xml_input)
        {

            try
            {
                xml_document input_file;
                xml_parse_result open_success = input_file.load_file(xml_input.c_str());
                if (!open_success) {
                    throw (int)errex::XML_READING_ERROR;
                }
                xml_node root = input_file.child("problem");


/*
    READ SIMULATION PARAMETERS.
*/
                mpreal dt;
                mpreal t;
                int precision_digits = 400;
                int output_digits = 15;
                const char_t* tmp;
                int AMin = -1;
                int AMax = -1;
                string output_location = "";

//..............Obtains the verbosity level from the input file.
                tmp = root.child("simulation_params").child("verbosity").child_value();
                tmp != "" ? __vbs__ = stoi(tmp) : __vbs__ = 0;

//..............Obtains the precision digits from the input file.
                tmp = root.child("simulation_params").child("precision_digits").child_value();
                tmp != "" ? precision_digits = stoi(tmp) : precision_digits = __dps__;

//..............Obtains the output precision digits from the input file.
                tmp = root.child("simulation_params").child("output_digits").child_value();
                tmp != "" ? output_digits = stoi(tmp) : output_digits = __dop__;


//..............IMPORTANT! Reads the precision before declarign any high-precision float vars.
                mpreal::set_default_prec(digits2bits(precision_digits));
                cout.precision(output_digits);

//..............Reads the substep interval in seconds.
                tmp = root.child("simulation_params").child("dt").child_value();
                tmp != "" ? dt = mpreal(tmp) : throw (int)errex::MISSING_SUBSTEP_SIZE;

//..............Reads the time step in seconds.
                tmp = root.child("simulation_params").child("time_step").child_value();
                tmp != "" ? t = mpreal(tmp) : throw (int)errex::MISSING_STEP_SIZE;

//..............Reads the final councentration output file location.
                tmp = root.child("simulation_params").child("output").child_value();
                tmp != "" ? output_location = tmp : output_location = ".//output.xml";
                if (__vbs__) cout << "Final nuclide concentrations will be written in " << output_location << endl;

//..............Reads the species names.
                if (!root.child("species")) throw (int)errex::MISSING_SPECIES_NAMES;
                vector<string> species_names;
                string species = root.child("species").child_value();
                if (strlen(root.child("species").attribute("amin").value()) > 0) {
                    xml_document source;
                    xml_parse_result open_success = source.load_file(root.child("species").attribute("source").value());
                    if (open_success)
                    {

                        if (strlen(root.child("species").attribute("amin").value()) > 0) {
                            AMin = stoi(root.child("species").attribute("amin").value());
                        }
                        else
                        {
                            cout << "warning <cnuctran::simulation::from_input()> AMin attribute is missing. Considering all nuclides with A > 0." << endl;
                            AMin = 0;
                        }
                        if (strlen(root.child("species").attribute("amax").value()) > 0) {
                            AMax = stoi(root.child("species").attribute("amax").value());
                        }
                        else
                        {
                            cout << "warning <cnuctran::simulation::from_input()> AMin attribute is missing. Considering all nuclides with A > 0." << endl;
                            AMax = 400;
                        }
                        species_names = get_nuclide_names(root.child("species").attribute("source").value(), AMin, AMax);
                        if (__vbs__) cout << "Building chains... A = [" << AMin << "," << AMax << "]. Total no. of nuclides = " << species_names.size() << endl;
                    }
                    else
                    {
                        throw (int)errex::MISSING_SPECIES_NAMES;
                    }
                }
                else
                {
                    stringstream ss = stringstream(species); string token;
                    while (getline(ss, token, ' '))
                    {
                        if (token != "")
                            species_names.push_back(token);
                    }
                    if (__vbs__) cout << "Building chains..." << endl;
                }


//..............Reads the initial concentrations.
                map<string, mpreal> w0;
                const char* w0_source = root.child("initial_concentration").attribute("source").value();
                const char* override_species_names = root.child("initial_concentration").attribute("override_species_names").value();
                if (w0_source != "")
                {
                    xml_document w0_doc;
                    xml_parse_result load_success = w0_doc.load_file(w0_source);
                    if (override_species_names == "true") species_names.clear();
                    if (load_success)
                    {
                        if (__vbs__) cout << "Reading the initial nuclide concentrations from " << w0_source << endl;
                        for (xml_node nuclide : w0_doc.child("nuclide_concentrations"))
                        {
                            string name = nuclide.attribute("name").value();
                            mpreal concentration = mpreal(nuclide.child_value());
                            if (override_species_names == "true")
                            {
                                species_names.push_back(name);
                                if (concentration != __zer__)
                                    w0[name] = concentration;
                            }
                            else
                            {
                                if (find(species_names.begin(), species_names.end(), name) != species_names.end())
                                    w0[name] = concentration;
                            }
                        }
                    }
                    else
                    {
                        throw (int)errex::MISSING_W0_SOURCE;
                    }
                }
                else
                {
                    
                    if (!root.child("initial_concentration").child_value()) throw (int)errex::MISSING_W0;
                    for (xml_node item : root.child("initial_concentration").children())
                    {
                        mpreal concentration = mpreal(item.child_value());
                        w0[item.attribute("species").value()] = concentration;
                        if (__vbs__ == 2) cout << "input<initial_concentration> species = " << item.attribute("species").value() << " w0 = " << concentration << endl;
                    }
                }
                


//..............Reads the rxn rates.
                map<string, map<string, mpreal>> rxn_rates;

                for (xml_node reaction : root.child("reaction_rates").children())
                {
                    mpreal rate = mpreal(reaction.child_value());
                    rxn_rates[reaction.attribute("species").value()][reaction.attribute("type").value()] = rate;
                    if (__vbs__ == 2) cout << "input<rxn_rates> species = " << reaction.attribute("species").value() << " type = " << reaction.attribute("type").value() << " rate = " << rate << endl;
                }


/*
    RUN THE SIMULATION.
*/
                solver sol = solver(species_names);
                build_chains(sol, rxn_rates, root.child("species").attribute("source").value());
                map<string, mpreal> w;
                w = sol.solve(w0, dt, t);

/*
    PRINTS THE OUTPUT.
*/
                // Prints to output file.
                stringstream ss("");
                ss << "<nuclide_concentrations amin=\"" << AMin << "\" amax=\"" << AMax << "\" n=\"" << sol.species_names.size() << "\" time_step=\"" << t << "\">" << endl;
                for (string species : sol.species_names)
                {
                    mpreal c = w[species];
                    if (c > __eps__)
                    {
                        ss << "\t" << "<nuclide name=\"" << species << "\">" << scientific << setprecision(output_digits) << w[species] << "</nuclide>" << endl;
                    }
                }
                ss << "</nuclide_concentrations>" << endl;
                ofstream file;
                file.open(output_location, ios::out);
                file << ss.str();
                file.close();
            }
/*
    EXCEPTIONS HANDLING.
*/
            catch (int e)
            {
                switch (e)
                {
                case errex::MISSING_SUBSTEP_SIZE:
                    cout << "fatal-error <cnuctran::simulation::from_input()> Substep size, dt, is not supplied." << endl;
                    exit(1);
                case errex::MISSING_STEP_SIZE:
                    cout << "fatal-error <cnuctran::simulation::from_input()> Time step, t, is not supplied." << endl;
                    exit(1);
                case errex::MISSING_SPECIES_NAMES:
                    cout << "fatal-error <cnuctran::simulation::from_input()> Species names are not supplied." << endl;
                    exit(1);
                case errex::MISSING_W0:
                    cout << "fatal-error <cnuctran::simulation::from_input()> Initial nuclide concentrations are not supplied." << endl;
                    exit(1);
                case errex::NUCLIDES_DATA_LOAD_FAILED:
                    cout << "fatal-error <cnuctran::simulation::from_input()> Could open source file." << endl;
                    exit(1);
                case errex::XML_READING_ERROR:
                    cout << "fatal-error <cnuctran::simulation::from_input()> An error has occurred while reading the XML input file: " << xml_input << endl;
                    exit(1);
                case errex::MISSING_W0_SOURCE:
                    cout << "fatal-error <cnuctran::simulation::from_input()> Could not open the initial nuclide concentrations XML file." << endl;
                    exit(1);
                    
                default:
                    cout << "fatal-error <cnuctran::simulation::from_input()> Unexpected error has occurred." << endl;
                    exit(1);
                }
            }


        }
    };
}





int main()
{

//..The program will only search for input.xml within the same directory.
    cnuctran::simulation::from_input(".\\input.xml");

    return 0;

}
