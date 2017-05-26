#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <assert.h>

using namespace std;

const int NS = 10;

double sign(double val) {
    return (0.0 < val) - (val < 0.0);
}

double calCst(const vector<double> &w, const vector<double> &v, const vector<int> &c, const double cl, const int i) {
    double cst = 0.5 * c[i] * pow(w[i]-v[i], 2);
    if (i == 0) {
        cst += cl * abs(w[i+1]-w[i]);
    } else if (i == v.size() - 1) {
        cst += cl * abs(w[i]-w[i-1]);
    } else {
        cst += cl * (abs(w[i]-w[i-1]) + abs(w[i+1]-w[i]));
    }
    return cst;
}

int main() {
    int n;
    double l;
    scanf("%d%lf", &n, &l);

    vector<double> ori;
    double d;
    for (int i = 0; i < n; i++) {
        scanf("%lf", &d);
        ori.emplace_back(d);
    }

    vector<double> v = ori;
    vector<double> w = ori;
    vector<int> c(n, 1);

    double osum = 0;
    for (int i = 0; i < n; i++) {
        if (i > 0) {
            osum += l * abs(w[i] - w[i-1]);
        }
    }

    for (int t = 0; t < NS; t++) {
        assert(v.size() == w.size() and w.size() == c.size());

        double cl = l / NS * (t+1);
        fprintf(stderr, "Iter: #%d with lambda %f\n", t, cl);

        double fsum = 0;
        for (int i = 0, k = 0; i < n; i++) {
            for (int j = 0; j < c[i]; j++, k++) {
                fsum += 0.5 * pow(w[i] - ori[k], 2);
            }
            if (i > 0) {
                fsum += cl * abs(w[i] - w[i-1]);
            }
        }
        
        bool update;
        do {
            update = false;

            for (size_t i = 0; i < v.size(); i++) {
                /* Step 1: Try descent cycle */
                if (v.size() >= 2) { // at least two element
                    double ocst = calCst(w, v, c, cl, i);
                    //fprintf(stderr, "ocst: %f (%f)\n", w[i], ocst);
                    if (i == 0) {
                        if ((c[i]*v[i]+cl)/c[i] <= w[i+1]) { // Case 1: w[i] <= w[i+1]
                            w[i] = (c[i]*v[i]+cl)/c[i];
                        } else if ((c[i]*v[i]-cl)/c[i] > w[i+1]) { // Case 2: w[i] > w[i+1]
                            w[i] = (c[i]*v[i]-cl)/c[i];
                        } else { // Case 3: try active constraint point: w[i+1]
                            w[i] = w[i+1];
                        }
                    } else if (i == v.size() - 1) {
                        if ((c[i]*v[i]+cl)/c[i] <= w[i-1]) { // Case 1: w[i] <= w[i-1]
                            w[i] = (c[i]*v[i]+cl)/c[i];
                        } else if ((c[i]*v[i]-cl)/c[i] > w[i-1]) { // Case 2: w[i] > w[i-1]
                            w[i] = (c[i]*v[i]-cl)/c[i];
                        } else { // Case 3: try active constraint point: w[i-1]
                            w[i] = w[i-1];
                        }
                    } else {
                        double wmax = max(w[i-1],w[i+1]);
                        double wmin = min(w[i-1],w[i+1]);
                        if ((c[i]*v[i]+2*cl)/c[i] <= wmin) { // Case 1: w[i] <= wmin
                            w[i] = (c[i]*v[i]+2*cl)/c[i];
                        } else if (v[i] > wmin and v[i] <= wmax) { // Case 2: wmin < w[i] <= wmax
                            w[i] = v[i];
                        } else if ((c[i]*v[i]-2*cl)/c[i] > wmax) { // Case 3: w[i] > wmax
                            w[i] = (c[i]*v[i]-2*cl)/c[i];
                        } else { // Case 4: try active constraint points: w[i-1], w[i+1]
                            // (1): 0.5 * c[i] * (wmin - v[i])^2 + |wmax-wmin| + |wmin-wmin|
                            // (2): 0.5 * c[i] * (wmax - v[i])^2 + |wmax-wmax| + |wmax-wmin|
                            if (pow(wmin-v[i],2) < pow(wmax-v[i],2)) {
                                w[i] = wmin; 
                            } else {
                                w[i] = wmax;
                            }
                        }
                    }
                        
                    double ncst = calCst(w, v, c, cl, i);
                    //fprintf(stderr, "ncst: %f (%f)\n", w[i], ncst);
                    if (ocst > ncst) { // cost reduce
                        //fprintf(stderr, "YESSS\n");
                        update = true;
                        continue; 
                    }
                }

                /* Step 2: Try fusion cycle */
                if (v.size() >= 3) { // at least three element
                    if (i == 0) {
                        // nothing can be done
                        continue;
                    }

                    double ocst = calCst(w, v, c, cl, i-1) + calCst(w, v, c, cl, i) - cl * abs(w[i]-w[i-1]);
                    //fprintf(stderr, "HELLO ocst: %f, %f (%f)\n", w[i-1], w[i], ocst);
                    double r; // w[i-1] = w[i] = r
                    if (i == 1) {
                        if ((c[i-1]*v[i-1]+c[i]*v[i]+cl)/(c[i-1]+c[i]) <= w[i+1]) { // Case 1: r <= w[i+1]
                            r = (c[i-1]*v[i-1]+c[i]*v[i]+cl)/(c[i-1]+c[i]);
                        } else if ((c[i-1]*v[i-1]+c[i]*v[i]-cl)/(c[i-1]+c[i]) > w[i+1]) { // Case 2: r > w[i+1]
                            r = (c[i-1]*v[i-1]+c[i]*v[i]-cl)/(c[i-1]+c[i]);
                        } else { // Case 3: try active constraint point: w[i+1]
                            r = w[i+1];
                        }
                    } else if (i == v.size() - 1) {
                        if ((c[i-1]*v[i-1]+c[i]*v[i]+cl)/(c[i-1]+c[i]) <= w[i-2]) { // Case 1: r <= w[i-2]
                            r = (c[i-1]*v[i-1]+c[i]*v[i]+cl)/(c[i-1]+c[i]);
                        } else if ((c[i-1]*v[i-1]+c[i]*v[i]-cl)/(c[i-1]+c[i]) > w[i-2]) { // Case 2: r > w[i-2]
                            r = (c[i-1]*v[i-1]+c[i]*v[i]-cl)/(c[i-1]+c[i]);
                        } else { // Case 3: try active constraint point: w[i-2]
                            r = w[i-2];
                        }
                    } else {
                        double wmax = max(w[i-2],w[i+1]);
                        double wmin = min(w[i-2],w[i+1]);
                        if ((c[i-1]*v[i-1]+c[i]*v[i]+2*cl)/(c[i-1]+c[i]) <= wmin) { // Case 1: r <= wmin
                            r = (c[i-1]*v[i-1]+c[i]*v[i]+2*cl)/(c[i-1]+c[i]); 
                        } else if ((c[i-1]*v[i-1]+c[i]*v[i])/(c[i-1]+c[i]) > wmin and 
                          (c[i-1]*v[i-1]+c[i]*v[i])/(c[i-1]+c[i]) <= wmax) { // Case 2: wmin < r <= wmax
                            r = (c[i-1]*v[i-1]+c[i]*v[i])/(c[i-1]+c[i]); 
                        } else if ((c[i-1]*v[i-1]+c[i]*v[i]-2*cl)/(c[i-1]+c[i]) > wmax) { // Case 3: r > wmax
                            r = (c[i-1]*v[i-1]+c[i]*v[i]-2*cl)/(c[i-1]+c[i]);
                        } else { // Case 4: try active constraint points: w[i-1], w[i+1]
                            // (1): 0.5 * c[i-1] * (wmin-v[i-1])^2 + 0.5 * c[i] * (wmin-v[i])^2 + |wmax-wmin| + |wmin-wmin|
                            // (2): 0.5 * c[i-1] * (wmax-v[i-1])^2 + 0.5 * c[i] * (wmax-v[i])^2 + |wmax-wmax| + |wmax-wmin|
                            if (c[i-1] * pow(wmin-v[i-1],2) + c[i] * pow(wmin-v[i],2) 
                              < c[i-1] * pow(wmax-v[i-1],2) + c[i] * pow(wmax-v[i],2)) {
                                r = wmin; 
                            } else {
                                r = wmax;
                            }
                        }
                    }

                    double wf = w[i-1];
                    double ws = w[i];
                    w[i-1] = w[i] = r;

                    double ncst = calCst(w, v, c, cl, i-1) + calCst(w, v, c, cl, i) - cl * abs(w[i]-w[i-1]);
                    //fprintf(stderr, "HELLO ncst: %f, %f (%f)\n", w[i-1], w[i], ncst);
                    if (ocst > ncst) { // cost reduce
                        update = true;
                        // keep updated w[i-1], w[i]
                        //fprintf(stderr, "HELLO YESSS\n");
                    } else {
                        // restore w[i-1], w[i]
                        //fprintf(stderr, "HELLO NOOOO\n");
                        w[i-1] = wf;
                        w[i] = ws;
                    }
                }
            }
        } while (update == true);

        /*
        fprintf(stderr, "WORLD\n");
        for (size_t i = 0; i < w.size(); i++) {
            fprintf(stderr, "%f ", w[i]);
        }
        fprintf(stderr, "\n");
        */

        /* Step 3: Try smoothing cycle */
        vector<double> nw;
        vector<double> nv;
        vector<int> nc;
        double last = w[0];
        double wsum = 0;
        double vsum = 0;
        double csum = 0;
        for (size_t i = 0; i < n; i++) {
            //fprintf(stderr, "%d %f %f %f\n", i, w[i], c[i], last);
            if (w[i] != last) {
                nw.emplace_back(wsum / csum);
                nv.emplace_back(vsum / csum);
                nc.emplace_back(csum);
                wsum = 0;
                vsum = 0;
                csum = 0;
                last = w[i];
            }
            wsum += c[i] * w[i];
            vsum += c[i] * v[i];
            csum += c[i];
        }
        nw.emplace_back(wsum / csum);
        nv.emplace_back(vsum / csum);
        nc.emplace_back(csum);

        w = nw;
        v = nv;
        c = nc;
        n = nw.size();

        double ssum = 0;
        double nsum = 0;
        for (int i = 0, k = 0; i < n; i++) {
            for (int j = 0; j < c[i]; j++, k++) {
                ssum += 0.5 * pow(w[i] - ori[k], 2);
                nsum += 0.5 * pow(w[i] - ori[k], 2);
            }
            if (i > 0) {
                ssum += cl * abs(w[i] - w[i-1]);
                nsum += l * abs(w[i] - w[i-1]);
            }
        }
        for (int i = 0; i < n; i++) {
            fprintf(stderr, "%f %f %d\n", w[i], v[i], c[i]);
        }
        fprintf(stderr, "Iter: #%d with reduce %f -> %f\n", t, fsum, ssum);
        fprintf(stderr, "Iter: #%d with reduce %f -> %f\n", t, osum, nsum);
    }
}
