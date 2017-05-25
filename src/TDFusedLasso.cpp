#include <cstdio>
#include <cstdlib>
#include <cfloat>
#include <cmath>
#include <vector>
#include <utility>
#include <set>
#include <algorithm>
#include "mex.h"

using namespace std;

const int NS = 10;

struct Node {
    Node(int r, int c, double v): loc(r, c), val(v) {}
    pair<int,int> loc;
    double val;
};

struct BigNode {
    BigNode(): v(0), w(0), edges(vector<pair<double,BigNode*>>()), nodes(vector<Node*>()) {}
    BigNode(Node* node): v(node->val), w(node->val) {
        nodes.emplace_back(node);
    }
    double v;
    double w;
    vector<pair<double,BigNode*>> edges;
    vector<Node*> nodes;
};

double calCst(const BigNode* bnode, const double w, const double cl) {
    double sum = 0;
    for (auto& p : bnode->edges) {
        sum += p.first * abs(w - p.second->w);
    }
    sum = 0.5 * bnode->nodes.size() * pow(w - bnode->v, 2) + cl * sum;
    return sum;
}

double calCst(const BigNode* bnode_f, const BigNode* bnode_s, const double w_f, const double w_s, const double cl) {
    double sum = 0;
    for (auto& p : bnode_f->edges) {
        //printf("s %f %f\n", p.first, abs(w_f - p.second->w));
        sum += p.first * abs(w_f - p.second->w);
    }
    for (auto& p : bnode_s->edges) {
        if (p.second == bnode_f) continue; // only add one-directional cost
        //printf("ss %f %f\n", p.first, abs(w_s - p.second->w));
        sum += p.first * abs(w_s - p.second->w);
    }
    /*
    printf("sum %f\n", cl * sum);
    printf("sum %f\n", 0.5 * bnode_f->nodes.size() * pow(w_f - bnode_f->v, 2));
    printf("sum %f\n", 0.5 * bnode_s->nodes.size() * pow(w_s - bnode_s->v, 2));
    */
    sum = cl * sum +
      0.5 * bnode_f->nodes.size() * pow(w_f - bnode_f->v, 2) +
      0.5 * bnode_s->nodes.size() * pow(w_s - bnode_s->v, 2);
    return sum;
}

double calCst(const vector<BigNode*>& vct, const double cl) {
    double sum = 0;
    for (size_t i = 0; i < vct.size(); i++) {
        for (size_t j = 0; j < vct[i]->nodes.size(); j++) {
            sum += 0.5 * pow(vct[i]->w - vct[i]->nodes[j]->val, 2);
        }
        for (auto& p : vct[i]->edges) {
            sum += 0.5 * cl * p.first * abs(vct[i]->w - p.second->w); // half of bi-directional cost
        }
    }
    return sum;
}

void dfs(BigNode* bnode, const double w, set<BigNode*>& vis, set<BigNode*>& svis) {
    vis.emplace(bnode);
    svis.emplace(bnode);
    for (auto& p : bnode->edges) {
        if (p.second->w == w and vis.find(p.second) == vis.end()) {
            dfs(p.second, w, vis, svis); 
        }
    }
}

void TDFusedLasso(double** v, int n, double lambda_1, double lambda_2, double** w) {
    vector<vector<BigNode*>> table(n, vector<BigNode*>(n, NULL));
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) { // upper triangle
            Node* node = new Node(i, j, v[i][j]); 
            table[i][j] = new BigNode(node);
        }
    }
    /* Add row-wise edges */
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n-1; j++) {
            table[i][j]->edges.emplace_back(1,table[i][j+1]);
            table[i][j+1]->edges.emplace_back(1,table[i][j]);
        }
    }
    /* Add col-wise edges */
    for (int j = 0; j < n; j++) {
        for (int i = 0; i <= j-1; i++) {
            table[i][j]->edges.emplace_back(1,table[i+1][j]);
            table[i+1][j]->edges.emplace_back(1,table[i][j]);
        }
    }

    vector<BigNode*> vct;
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            vct.emplace_back(table[i][j]);
        }
    }
    
    for (size_t i = 0; i < vct.size(); i++) {
        //mexPrintf("%f %f\n", vct[i]->w, vct[i]->v);
    }

    double osum = calCst(vct, lambda_2);

    /* Increase lambda slightly */
    for (int t = 0; t < NS; t++) {
        double cl = lambda_2 / NS * (t+1); 
        //mexPrintf("Iter: #%d with lambda %f\n", t, cl);

        double fsum = calCst(vct, cl);

        bool update;
        do {
            update = false;

            for (size_t i = 0; i < vct.size(); i++) {
                /* Step 1: Try descent cycle */
                {
                    double ocst = calCst(vct[i], vct[i]->w, cl);
                    //mexPrintf("step1: #%d ocst %f\n", int(i), ocst);

                    vector<double> wlst;
                    for (auto& p : vct[i]->edges) {
                        wlst.emplace_back(p.second->w);
                    }
                    sort(wlst.begin(), wlst.end());
                    wlst.emplace_back(DBL_MAX); // add last boundary

                    bool descent_flag = false;
                    /* Case 1: piecewise linear solution */
                    double nw;
                    for (size_t j = 0; j < wlst.size(); j++) {
                        nw = 0;
                        for (auto& p : vct[i]->edges) {
                            if (p.second->w >= wlst[j]) {
                                nw += p.first;
                            } else {
                                nw -= p.first;
                            }
                        }
                        nw = (vct[i]->v * vct[i]->nodes.size() + cl * nw) / vct[i]->nodes.size();

                        if ((j == 0 and nw < wlst[j]) or (wlst[j-1] <= nw and nw < wlst[j])) {
                            descent_flag = true;
                            break;
                        }
                    }

                    /* Case 2: active constraint points */
                    if (descent_flag == false) {
                        double minc = DBL_MAX;
                        double minw = -1;
                        wlst.pop_back(); // remove last boundary
                        for (auto& w : wlst) {
                            double cst = calCst(vct[i], w, cl);
                            if (cst < minc) {
                                minc = cst;
                                minw = w;
                            }
                        }
                        nw = minw;
                    }

                    /* Joint consideration for Case 1 and Case 2 */
                    double ncst = calCst(vct[i], nw, cl);
                    //mexPrintf("step1: #%d ncst %f\n", int(i), ncst);
                    if (ocst > ncst) { // cost reduce
                        update = true;
                        vct[i]->w = nw;
                        continue;
                    } else {
                    }
                }

                /* Step 2: Try fusion cycle */
                {
                    for (auto& p : vct[i]->edges) {
                        double ocst = calCst(vct[i], p.second, vct[i]->w, p.second->w, cl);
                        //mexPrintf("step2: #%d ocst %f\n", int(i), ocst);

                        vector<pair<double,BigNode*>> medges;
                        for (auto& q : vct[i]->edges) {
                            if (q.second == p.second) continue;
                            medges.emplace_back(q);
                        }
                        for (auto& q : p.second->edges) {
                            if (q.second == vct[i]) continue;
                            medges.emplace_back(q);
                        }

                        vector<double> wlst;
                        for (auto& q: medges) {
                            wlst.emplace_back(q.second->w);
                        }
                        sort(wlst.begin(), wlst.end());
                        wlst.emplace_back(DBL_MAX); // add last boundary

                        bool descent_flag = false;
                        /* Case 1: piecewise linear solution */
                        double nw;
                        for (size_t j = 0; j < wlst.size(); j++) {
                            nw = 0;
                            for (auto& q : medges) {
                                if (q.second->w >= wlst[j]) {
                                    nw += q.first;
                                } else {
                                    nw -= q.first;
                                }
                            }
                            nw = (vct[i]->v * vct[i]->nodes.size() + p.second->v * p.second->nodes.size() + cl * nw) 
                              / (vct[i]->nodes.size() + p.second->nodes.size());

                            if ((j == 0 and nw < wlst[j]) or (wlst[j-1] <= nw and nw < wlst[j])) {
                                descent_flag = true;
                                break;
                            }
                        }

                        /* Case 2: active constraint points */
                        if (descent_flag == false) {
                            double minc = DBL_MAX;
                            double minw = -1;
                            wlst.pop_back(); // remove last boundary
                            for (auto& w : wlst) {
                                double cst = calCst(vct[i], p.second, w, w, cl);
                                if (cst < minc) {
                                    minc = cst;
                                    minw = w;
                                }
                            }
                            nw = minw;
                        }

                        /* Joint consideration for Case 1 and Case 2 */
                        double ncst = calCst(vct[i], p.second, nw, nw, cl);
                        //mexPrintf("step2: #%d ncst %f\n", int(i), ncst);
                        if (ocst > ncst) { // cost reduce
                            update = true;
                            //mexPrintf("HI!!!!!\n");
                            vct[i]->w = nw;
                            p.second->w = nw;
                        }
                    }
                }
            }
        } while (update);

        /* Step 3: Try smoothing cycle */
        vector<BigNode*> nvct;
        set<BigNode*> vis;
        for (size_t i = 0; i < vct.size(); i++) {
            if (vis.find(vct[i]) != vis.end()) continue;

            /* Traverse whole graph with same w */
            set<BigNode*> svis;
            dfs(vct[i], vct[i]->w, vis, svis);

            /* Merge groups */
            BigNode* nnode = new BigNode();
            for (auto& bnode : svis) {
                nnode->w += bnode->w * bnode->nodes.size();
                nnode->v += bnode->v * bnode->nodes.size();
                for (auto& node : bnode->nodes) {
                    nnode->nodes.emplace_back(node);
                }
            }
            nnode->w /= nnode->nodes.size();
            nnode->v /= nnode->nodes.size();

            nvct.emplace_back(nnode);
        }

        /* Build new edges */
        for (size_t i = 0; i < nvct.size(); i++) {
            for (size_t j = i+1; j < nvct.size(); j++) {
                int cnt = 0;
                for (auto& nu : nvct[i]->nodes) {
                    for (auto& nv : nvct[j]->nodes) {
                        if (abs(nu->loc.first-nv->loc.first) + abs(nu->loc.second-nv->loc.second) == 1) {
                            cnt += 1;
                        }
                    }
                }

                if (cnt > 0) {
                    nvct[i]->edges.emplace_back(cnt, nvct[j]);
                    nvct[j]->edges.emplace_back(cnt, nvct[i]);
                }
            }
        }

        /* Delete old objects */
        for (size_t i = 0; i < vct.size(); i++) {
            delete vct[i];
        }
        vct = nvct;

        for (size_t i = 0; i < vct.size(); i++) {
            //mexPrintf("%f %f %d\n", vct[i]->w, vct[i]->v, int(vct[i]->nodes.size()));
            for (auto& node : vct[i]->nodes) {
                //mexPrintf("%f ", node->val); 
            }
            //mexPrintf("\n");
        }

        double ssum = calCst(vct, cl);
        double nsum = calCst(vct, lambda_2);
        //mexPrintf("Iter: #%d with reduce %f -> %f\n", t, fsum, ssum);
        //mexPrintf("Iter: #%d with reduce %f -> %f\n", t, osum, nsum);
    }   

    for (size_t i = 0; i < vct.size(); i++) {
        for (auto& node : vct[i]->nodes) {
            double res = vct[i]->w;
            if (res > lambda_1) {
                res -= lambda_1;
            } else if (res < -lambda_1) {
                res += lambda_1;
            } else {
                res = 0;
            }
            w[node->loc.first][node->loc.second] = res;
        }
    }

    for (size_t i = 0; i < vct.size(); i++) {
        for (auto& node : vct[i]->nodes) {
            delete node;
        }
        delete vct[i];
    }
}

void mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    /* Input arguments */
    double* v = mxGetPr(prhs[0]);
    int n = (int)mxGetScalar(prhs[1]);
    double lambda_1 = mxGetScalar(prhs[2]);
    double lambda_2 = mxGetScalar(prhs[3]);
    
    double** fv = new double*[n];
    for (int i = 0; i < n; i++) {
        fv[i] = new double[n]();
    }
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            fv[i][j] = v[j*n + i];
        }
    }

    double** fw = new double*[n];
    for (int i = 0; i < n; i++) {
        fw[i] = new double[n]();
    }
    TDFusedLasso(fv, n, lambda_1, lambda_2, fw);

    /* Output arguments */
    plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
    double* w = mxGetPr(plhs[0]);
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            w[j*n + i] = fw[i][j];
        }
    }

    /* Release memory */
    for (int i = 0; i < n; i++) {
        delete fv[i];
        delete fw[i];
    }
    delete[] fv;
    delete[] fw;
}
