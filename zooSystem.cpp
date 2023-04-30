//Mariia Shmakova
//m.shmakova@innopolis.university

#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace std;

auto v(int v0, int k0, double alpha1, double alpha2, double beta1, double beta2, double time){
    return(v0* cos(sqrt(alpha1*alpha2)*time) - k0*((sqrt(alpha2)*beta1)/(beta2*(sqrt(alpha1))))
    * sin(sqrt(alpha1*alpha2)*time));
}

auto k(int v0, int k0, double alpha1, double alpha2, double beta1, double beta2, double time){
    return(v0* ((sqrt(alpha1)*beta2)/(beta1*(sqrt(alpha2))))*sin(sqrt(alpha1*alpha2)*time)
    + k0*cos(sqrt(alpha1*alpha2)*time));
}


int main(){

    int v0;
    int k0;
    double alpha1;
    double beta1;
    double alpha2;
    double beta2;
    double time;
    int n;
    double variable=0.0;

    cin >> v0;
    cin >> k0;
    cin >> alpha1;
    cin >> beta1;
    cin >> alpha2;
    cin >> beta2;
    cin >> time;
    cin >> n;

    double voOld= v0;
    double koOld = k0;

    v0 -= alpha2/beta2;
    k0 -= alpha1/beta1;

    cout << "t:\n";

    while(variable <= time){
        cout << fixed << setprecision(2);
        if(variable == time){
            cout << variable;
        }else{
            cout << variable << " ";
        }

        variable += time / n;
    }


    variable = 0.0;
    cout << fixed << setprecision(2);
    cout << "\nv:\n" << voOld << " ";
    variable += time / n;

    while(variable <= time){
        cout << fixed << setprecision(2);
        if(variable == time){
            cout << v(v0, k0, alpha1, alpha2, beta1, beta2, variable) + alpha2/beta2;
        }else{
            cout << v(v0, k0, alpha1, alpha2, beta1, beta2, variable) + alpha2/beta2 << " ";
        }

        variable += time / n;
    }


    variable = 0.0;
    cout << fixed << setprecision(2);
    cout << "\nk:\n" << koOld << " ";
    variable += time / n;

    while(variable <= time){
        cout << fixed << setprecision(2);
        if(variable == time){
            cout << k(v0, k0, alpha1, alpha2, beta1, beta2, variable) + alpha1/beta1;
        }else{
            cout << k(v0, k0, alpha1, alpha2, beta1, beta2, variable) + alpha1/beta1 << " ";
        }

        variable += time / n;
    }
    return 0;

}
