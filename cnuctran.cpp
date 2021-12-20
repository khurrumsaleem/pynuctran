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


    const mpreal __two__ = mpreal("2.0");
    const mpreal __one__ = mpreal("1.0");
    const mpreal __neg__ = mpreal("-1.0");
    const mpreal __zer__ = mpreal("0.0");
    const mpreal __eps__ = mpreal("1e-30");
    bool verbosity = true;
    const int    __nop__ = -1;

    void row_operation(int& irow,
        map<int, map<int, mpreal>>& sd,
        map<int, map<int, mpreal>>& od,
        map<int, map<int, mpreal>>& rd)
    {
        map<int, mpreal>* rdr = &rd[irow];
        map<int, mpreal>* sdd = &sd[irow];

        for (map<int, mpreal>::iterator it1 = sdd->begin(); it1 != sdd->end(); it1++) {
            int icol = it1->first;
            map<int, mpreal>* odd = &od[icol];
            mpreal x = (*sdd)[icol];
            for (map<int, mpreal>::iterator it2 = odd->begin(); it2 != odd->end(); it2++) {
                int ocol = it2->first;
                (*rdr)[ocol] += x * (*odd)[ocol];
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

        smatrix multiply(smatrix& other)
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
            r = this->multiply(other);
            return r;
        }

        smatrix pow(mpz_t n)
        {

            smatrix result = this->copy();
            mpz_t one; mpz_init(one); mpz_set_ui(one, 1);
            mpz_t two; mpz_init(two); mpz_set_ui(two, 1);

            if (mpz_cmp(n, one) == 0)
                return result;

            mpz_t int_div; mpz_init(int_div); mpz_fdiv_q(int_div, n, two);
            result = result.pow(int_div);

            mpz_t mod_res; mpz_init(mod_res); mpz_mod(mod_res, n, two);
            if (mod_res)
                return result * result * *this;
            else
                return result * result;

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
                this->G.push_back(vector<vector<int>>({ {__nop__} }));
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

            cout << "define removal: parent = " << this->species_names[species_index];
            cout << "\tindex = " << species_index;
            cout << "\trate = " << rate;
            cout << "\tproducts = ";


            int i = 0;
            for (int p : products)
            {
                if (p != __nop__)
                    cout << "[" << p << "]\t" << this->species_names[p];
                else
                {
                    cout << "[" << p << "]\tnot-tracked";
                    continue;
                }

                if (fission_yields.size() != 0)
                    cout << " (fission_yield = " << fission_yields[i++] << ")\t";
                else
                    cout << "\t";
            }
            cout << endl;


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
            vector<vector<mpreal>> A = vector<vector<mpreal>>();
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
                        mpreal kron = mpreal(to_string((int)(l == i)).c_str());
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
            mpreal t,
            mpz_t substeps)
        {
            vector<vector<mpreal>> w0_matrix = vector<vector<mpreal>>();
            for (int i = 0; i < this->__I__; i++)
                w0_matrix.push_back(vector<mpreal>({ __zer__ }));
            for (int i = 0; i < this->__I__; i++)
                if (w0.count(this->species_names[i]) == 1)
                    w0_matrix[i][0] = w0[this->species_names[i]];
            smatrix converted_w0 = smatrix(w0_matrix);

            mpreal dt = t / substeps;
            smatrix A = this->prepare_transfer_matrix(dt);
            cout << "Done building transfer matrix. Size = (" << A.shape.first << ", " << A.shape.second << ")" << endl;

            smatrix An = A.pow(substeps);
            cout << "Done computing sparse matrix power." << endl;

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
                    decay_rate = log(__two__) / mpreal(species.attribute("half_life").value());
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
                                            stringstream ss(data.value()); string token;
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
                                                        stringstream ss(param.value()); string token;
                                                        vector<string> products;
                                                        while (getline(ss, token, ' '))
                                                            if (token != "")
                                                                products.push_back(token);
                                                    }

                                                    if (string(param.name()) == "data")
                                                    {
                                                        stringstream ss(param.value()); string token;
                                                        while (getline(ss, token, ' '))
                                                            yields.push_back(mpreal(token));
                                                    }



                                                }

                                                mpreal total_fission_rate = rxn_rates[parent]["fission"];
                                                vector<mpreal> yields_to_add;
                                                vector<int> daughters_id_to_add;

                                                for (string product : products)
                                                {
                                                    vector<string>::iterator it_product;
                                                    it_product = std::find(species_names.begin(), species_names.end(), product);
                                                    if (it_product != species_names.end())
                                                    {
                                                        int product_id = distance(species_names.begin(), it_product);
                                                        daughters_id_to_add.push_back(product_id);
                                                        it_product = std::find(products.begin(), products.end(), product);
                                                        product_id = distance(products.begin(), it_product);
                                                        yields_to_add.push_back(yields[product_id]);
                                                    }
                                                    it_parent = std::find(species_names.begin(), species_names.end(), parent);
                                                    parent_id = distance(products.begin(), it_parent);
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
    };
}

int main()
{
    using namespace cnuctran;
    mpreal::set_default_prec(digits2bits(600));
    cout.precision(10);

    solver sol = solver(depletion_scheme::get_nuclide_names("E:\\chain_endfb71.xml", 200,300));
    map<string, map<string, mpreal>> rxn_rates;
    depletion_scheme::build_chains(sol, rxn_rates, "E:\\chain_endfb71.xml");


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
