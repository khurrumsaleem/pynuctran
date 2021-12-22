// cnuctran.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <iomanip>
#include <mpfr.h>
#include <mpir.h>
#include <mpreal.h>
#include <map>
#include <vector>
#include "pugixml.hpp"

using namespace mpfr;
using namespace std;
using namespace pugi;

namespace cnuctran {


   
    const mpreal __two__ = mpreal("2.0", digits2bits(500));
    const mpreal __one__ = mpreal("1.0", digits2bits(500));
    const mpreal __neg__ = mpreal("-1.0", digits2bits(500));
    const mpreal __zer__ = mpreal("0.0", digits2bits(500));
    const mpreal __eps__ = mpreal("1e-30", digits2bits(500));
    bool verbosity = true;
    const int    __nop__ = -1;

    map<int, mpreal>* rdr;
    map<int, mpreal>* sdd;
    map<int, mpreal>* odd;

    void row_operation(int& irow,
        map<int, map<int, mpreal>>& sd,
        map<int, map<int, mpreal>>& od,
        map<int, map<int, mpreal>>& rd)
    {
        rdr = &rd[irow];
        sdd = &sd[irow];

        for (map<int, mpreal>::iterator it1 = sdd->begin(); it1 != sdd->end(); it1++) {
            int icol = it1->first;
            odd = &od[icol];
            mpreal x = (*sdd)[icol];
            for (map<int, mpreal>::iterator it2 = odd->begin(); it2 != odd->end(); it2++) {
                int ocol = it2->first;
                (*rdr)[ocol] += x * (*odd)[ocol];
                if ((*rdr)[ocol] < __eps__)
                    (*rdr)[ocol] = __zer__;
            }
        }

        return;

    }

    class smatrix {

    public:


        pair<int, int> shape;
        map<int, map<int, mpreal>> data;

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
            map<int, map<int, mpreal>>* sd = &this->data;
            map<int, map<int, mpreal>>* od = &other.data;
            map<int, map<int, mpreal>>* rd = &result.data;

            for (int row = 0; row < sx; row++) {
                row_operation(row, *sd, *od, *rd);
            }

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

           /*
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
           
           */
       
            
            



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
                    cout << "FATAL-ERROR <cnuctran::solver::add_removal(...)> Insufficient fission yields given for species " <<
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
                cout << "FATAL-ERROR <cnuctran::solver::add_removal(...)> Invalid removal definition for isotope " <<
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
            mpz_t substeps)
        {
            vector<vector<mpreal>> w0_matrix;
            for (int i = 0; i < this->__I__; i++)
                w0_matrix.push_back(vector<mpreal>({ __zer__ }));
            for (int i = 0; i < this->__I__; i++)
                if (w0.count(this->species_names[i]) == 1)
                    w0_matrix[i][0] = w0[this->species_names[i]];
            smatrix converted_w0 = smatrix(w0_matrix);

            cout << "The substep interval was set to " << dt << " secs. " << endl;
            smatrix A = this->prepare_transfer_matrix(dt);
            cout << "Done building transfer matrix. Size = (" << A.shape.first << ", " << A.shape.second << ")" << endl;
            smatrix An = A.pow(substeps);
            cout << "Done computing sparse matrix power." << endl;
            //cout << A.to_string(100) << endl;
            smatrix w = An * converted_w0;
            cout << "Done computing concentrations." << endl;

            map<string, mpreal> out;
            for (int i = 0; i < this->__I__; i++)
                out[this->species_names[i]] = w.data[i][0];

            return out;

        }
    };

    class depletion_scheme
    {

    public:

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
            cout << "Done building chains." << endl;
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

        static void simulate_from_input(string xml_input)
        {
            xml_document input_file;
            xml_parse_result open_success = input_file.load_file(xml_input.c_str());
            if (!open_success) {
                cout << "FATAL-ERROR <cnuctran::depletion_scheme::read_input(...)> An error has occurred while reading the XML input file." << endl;
                exit(1);
            }
            xml_node root = input_file.child("problem");



            // Read simulation parameters.
            mpreal dt;
            mpz_t substeps;
            int precision_digits = 400;
            int output_digits = 15;


            precision_digits = stoi(root.child("simulation_params").child("precision_digits").child_value());
            cout << "input<simulation_params> precision_digits = " << precision_digits << endl;

            output_digits = stoi(root.child("simulation_params").child("output_digits").child_value());
            cout << "input<simulation_params> output_digits = " << output_digits << endl;

            // Set the precision...
            mpreal::set_default_prec(digits2bits(precision_digits));
            cout.precision(output_digits);

            // Set dt and substeps...
            dt = mpreal(root.child("simulation_params").child("dt").child_value());
            cout << "input<simulation_params> dt = " << dt << endl;

            mpz_init(substeps);
            mpz_set_str(substeps, root.child("simulation_params").child("substeps").child_value(), 10);
            //substeps = mpreal(root.child("simulation_params").child("substeps").child_value());
            cout << "input<simulation_params> substeps = " << substeps << endl;

            //Read the nuclides.
            vector<string> species_names;
            string species = root.child("species").child_value();
            if (strlen(root.child("species").attribute("AMin").value()) > 0) {
                xml_document source;
                xml_parse_result open_success = source.load_file(root.child("species").attribute("source").value());
                if (open_success)
                {
                    int AMin = -1;
                    int AMax = -1;
                    if (strlen(root.child("species").attribute("AMin").value()) > 0) {
                        AMin = stoi(root.child("species").attribute("AMin").value());
                    }
                    if (strlen(root.child("species").attribute("AMax").value()) > 0) {
                        AMax = stoi(root.child("species").attribute("AMax").value());
                    }
                    species_names = depletion_scheme::get_nuclide_names(root.child("species").attribute("source").value(), AMin, AMax);
                    cout << "Done loading species names, A = [" << AMin << "," << AMax << "]" << endl;
                }
                else
                {
                    cout << "fatal-error <cnuctran::depletion_scheme::simulate_from_input(...)> Cannot open source file." << endl;
                    exit(1);
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
                cout << "Done loading species names." << endl;
            }
            

            //Read the initial concentration.
            map<string, mpreal> w0;
            for (xml_node item : root.child("initial_concentration").children())
            {
                mpreal concentration = mpreal(item.child_value());
                w0[item.attribute("species").value()] = concentration;
                cout << "input<initial_concentration> species = " << item.attribute("species").value() << " w0 = " << concentration << endl;
            }
        

            //Read the rxn rates.
            map<string, map<string, mpreal>> rxn_rates;

            for (xml_node reaction : root.child("reaction_rates").children())
            {
                mpreal rate = mpreal(reaction.child_value());
                rxn_rates[reaction.attribute("species").value()][reaction.attribute("type").value()] = rate;
                cout << "input<rxn_rates> species = " << reaction.attribute("species").value() << " type = " << reaction.attribute("type").value() << " rate = " << rate << endl;
            }



            solver sol = solver(species_names);
            depletion_scheme::build_chains(sol, rxn_rates, root.child("species").attribute("source").value());
            map<string, mpreal> w = sol.solve(w0, dt, substeps);

            for (string species : sol.species_names)
            {
                mpreal c = w[species];
                if (c > __eps__)
                    cout << species << "\t" << w[species] << endl;
            }
            
        }
    };
}

int main()
{
    using namespace cnuctran;

    depletion_scheme::simulate_from_input("E:\\input.xml");

    exit(0);
    


}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
