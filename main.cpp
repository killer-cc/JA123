#include <iostream>
#include <cmath>
#include <random>
#include <cstdlib>
#include <map>
#include <vector>
#include <iomanip>
#include <fstream>

using namespace std;
fstream f ("123.epin",ios::out);
int epin_count=1;
// Use a function pointer instead of constant function for extendability
// Example: std::function<double(float)> FitnessFunction;
double function(float a) {
    return 20.0+double(exp(1))-20.0*double(exp(-0.2*double(sqrt(double(pow(a,2.0))))))-double(exp(double(cos(2.0*double(acos(-1.0))*a))));
    //return 10.0 + (pow(a, 2.0) - (10.0 * double(cos(2.0 * double(acos(-1.0)) * a))));
    //return fabs(a);
    //return a*a;
    //return fabs(a+5.281);
    /*if (a <= 0) {
        return pow(10, -a);
    } else if (a > 0 && a <= 20) {
        return pow(0.999, a);
    } else {
        return pow(2, a);
    }*/
}
float best(float a, float b, float c){
    float best;
    best=function(a)<function(b)?a:b;
    best=function(best)<function(c)?best:c;
    return best;
}
float min(vector<float> a){
    float min=a.at(0);
    for(int i=0;i<a.size();i++){
        if(function(a.at(i))<function(min))
            min=a.at(i);
    }
    return min;
}
class Jaguar {
public:
    /*
     vector<float> territory;
    const float up_bound = 100;//上界
    const float low_bound = -100;//下界
    float punish = 2147483647.0;//很差的值
    float place;//當前位置
    double fitness;//當前值
    float past_best_place;//過去最好位置
    double past_best_fitness;//過去最好位置的值
    float step;//狩獵步伐
    bool right_or_left;//狩獵左右
    bool hunt_check ;//判斷狩獵結束的check，結束為true
    int acc_step_count;//加速數值
    float access ;//總存取數
    float jumping_radius;//跳躍半徑
    float now_territory;//現在所在坑
     */
    vector<float> territory;
    const float up_bound = 32.768;
    const float low_bound = -32.768;
    float punish = 2147483647.0;
    float place;
    double fitness;
    float past_best_place;
    double past_best_fitness;
    float step;
    bool right_or_left;
    bool hunt_check ;
    int acc_step_count;
    float access ;
    float jumping_radius;
    float now_territory;
    float epin_radius;
    int final_access;

    Jaguar();

    void init();

    void init(float);

    void look_around();

    void hunting_jump();

    float speed_down_around();

    void hunting();

    float hunting(float);

    void jumping();
};
Jaguar::Jaguar() {
    place = 0;
    fitness = 0;
}
void Jaguar::init() {
    right_or_left = false;
    hunt_check = false;
    acc_step_count = 1;
    jumping_radius=0;
    now_territory=0;
    place = rand() / float(RAND_MAX) * (up_bound - low_bound) + low_bound;
    //cout << setprecision(20) << place << endl;
    jumping_radius=place;
    past_best_place = place;
    fitness = function(place);
    f<<setprecision(20)<<"*"<<epin_count<<" "<<fitness<<": "<<place<<",10,0,12"<<endl;
    epin_count++;
    past_best_fitness = fitness;
    access++;
    step = powf(2.0, int(log2(up_bound)) - 11);
}
void Jaguar::init(float a) {
    right_or_left = false;
    hunt_check = false;
    acc_step_count = 1;
    jumping_radius=0;
    place = a;
    //cout << setprecision(20) << place << endl;
    past_best_place = place;
    jumping_radius=place;
    fitness = function(place);
    f<<setprecision(20)<<"*"<<epin_count<<" "<<fitness<<": "<<place<<",10,"<<epin_radius<<",12"<<endl;
    epin_count++;
    past_best_fitness = fitness;
    access++;
    step = powf(2.0, int(log2(up_bound)) - 11);
}
void Jaguar::look_around() {
    float place_l, place_r;
    double fitness_l, fitness_r;
    place_l = place - step * acc_step_count;
    place_r = place + step * acc_step_count;
    if (place_r < low_bound || place_r > up_bound)
        fitness_r = punish;
    else{
        fitness_r = function(place_r);
        access++;
    }
    //cout << setprecision(20) << place_r << "," << step * acc_step_count<< "," << access << ","<< fitness_r << endl;
    f<<setprecision(20)<<"*"<<epin_count<<" "<<fitness_r<<": "<<place_r<<",10,0,12"<<endl;
    epin_count++;
    if (place_l < low_bound || place_l > up_bound)
        fitness_l = punish;
    else{
        fitness_l = function(place_l);
        access++;
    }
    //cout << setprecision(20) << place_l << "," << step * acc_step_count << "," << access << ","<< fitness_l << endl;
    f<<setprecision(20)<<"*"<<epin_count<<" "<<fitness_l<<": "<<place_l<<",10,0,12"<<endl;
    epin_count++;
    place=fitness<fitness_l?place:place_l;
    fitness = fitness<fitness_l?fitness:fitness_l;
    place = fitness < fitness_r ? place : place_r;
    fitness = fitness < fitness_r ? fitness : fitness_r;
    if (past_best_fitness > fitness) {
        past_best_place = place;
        past_best_fitness = fitness;
    }
    right_or_left = fitness_l < fitness_r;
}
void Jaguar::hunting_jump() {
    bool stop_flag = false;
    int dir=0;
    if(right_or_left)
        dir=-1;
    else
        dir=1;
    while (!stop_flag) {
        acc_step_count = acc_step_count * 2;
        float vr_place;
        double vr_fitness;
        vr_place = place + step * acc_step_count*dir;
        if (vr_place > low_bound && vr_place < up_bound) {
            vr_fitness = function(vr_place);
            access++;
        } else {
            vr_fitness = punish;
        }
        //cout << setprecision(20) << vr_place << "," << step * acc_step_count << "," << access << ","<< vr_fitness << endl;
        f<<setprecision(20)<<"*"<<epin_count<<" "<<vr_fitness<<": "<<vr_place<<",10,0,12"<<endl;
        epin_count++;
        if (vr_fitness < fitness) {
            place = vr_place;
            fitness = vr_fitness;
            past_best_place = place;
            past_best_fitness = fitness;
        } else {
            acc_step_count = acc_step_count / 2;
            vr_place = place + step * acc_step_count*dir;
            if (vr_place > low_bound && vr_place < up_bound) {
                vr_fitness = function(vr_place);
                access++;
            } else {
                vr_fitness = punish;
            }
            place = vr_place;
            fitness = vr_fitness;
            if (past_best_fitness > fitness) {
                past_best_place = place;
                past_best_fitness = fitness;
            }
            //cout << setprecision(20) << place << "," << step * acc_step_count << "," << access<< "," << fitness << endl;
            f<<setprecision(20)<<"*"<<epin_count<<" "<<fitness<<": "<<place<<",10,0,12"<<endl;
            epin_count++;
            stop_flag = true;
        }
    }
}
float Jaguar::speed_down_around() {
    while (acc_step_count > 1) {
        acc_step_count = acc_step_count / 2;
        place = past_best_place;
        fitness = past_best_fitness;
        look_around();
    }
    double a = past_best_fitness;
    while (past_best_fitness == a) {
        place = past_best_place;
        fitness = past_best_fitness;
        step = step / 2;
        if ((place+step)==(place-step)) {
            jumping_radius=fabs(jumping_radius-place);
            hunt_check=true;
            return place;
        }
        else {
            look_around();
        }
    }
    return 1;
}
void Jaguar::hunting() {
    init();
    float a=place;
    look_around();
    float final_place=1;
    while(a==place){
        final_place=speed_down_around();
        if(hunt_check){
            break;
        }
    }
    while (!hunt_check) {
        hunting_jump();
        final_place=speed_down_around();
    }
    if(!(find(territory.begin(),territory.end(),place)!=territory.end())){
        territory.push_back(final_place);
    }
}
float Jaguar::hunting(float a) {
    hunt_check= false;
    float final_place=1;
    init(a);
    float b=place;
    look_around();
    while(b==place){
        final_place=speed_down_around();
        if(hunt_check){
            break;
        }
    }
    while (!hunt_check) {
        hunting_jump();
        final_place=speed_down_around();
    }
    return final_place;
}
void Jaguar::jumping(){
    now_territory=territory.at(0);
    final_access=access;
    float territory_r=2147483647.0,territory_l=2147483647.0;
    int count=0;
    float jumping_place_r,jumping_place_l;
    float this_time_radius=jumping_radius;
    bool outside_check_r= false,outside_check_l= false;
    bool new_check_r=false,new_check_l=false;
    float sw_this_time_radius_r=this_time_radius*2.0f,sw_this_time_radius_l=this_time_radius*2.0f;
    do{

        if(count!=0){
            float a;
            //cout<<fabs(jumping_place_r-now_territory)<<endl;
            //cout<<fabs(jumping_place_l-now_territory)<<endl;
            if(!outside_check_r && !new_check_r&&!outside_check_l && !new_check_l&&fabs(jumping_place_r-now_territory)>=fabs(jumping_place_l-now_territory)){
                a=fabs(jumping_place_r-now_territory);
            }
            else{
                a=fabs(jumping_place_l-now_territory);
            }
            if(outside_check_r || new_check_r)
            {
                a=fabs(jumping_place_l-now_territory);
            }
            if(outside_check_l || new_check_l){
                a=fabs(jumping_place_r-now_territory);
            }
            sw_this_time_radius_r=a*2.0f;
            sw_this_time_radius_l=a*2.0f;
        }
        epin_radius=sw_this_time_radius_r;
        jumping_place_r=territory.at(0)+sw_this_time_radius_r;
        if(!outside_check_r && jumping_place_r>=up_bound){
            territory_r=hunting(up_bound);
            if(!(find(territory.begin(),territory.end(),territory_r)!=territory.end())){
                territory.push_back(territory_r);
                final_access=access;
                new_check_r=true;
            }
            outside_check_r= true;
        }
        if (!outside_check_r && !new_check_r){
            territory_r=hunting(jumping_place_r);
            if(!(find(territory.begin(),territory.end(),territory_r)!=territory.end())){
                territory.push_back(territory_r);
                final_access=access;
                new_check_r=true;
            }
        }
        //cout<<territory_r<<endl;
        epin_radius=sw_this_time_radius_l;
        jumping_place_l=territory.at(0)-sw_this_time_radius_l;
        if(!outside_check_l && jumping_place_l<=low_bound){
            territory_l=hunting(low_bound);
            if(!(find(territory.begin(),territory.end(),territory_l)!=territory.end())){
                territory.push_back(territory_l);
                final_access=access;
                new_check_l=true;
            }
            outside_check_l= true;
        }
        if (!outside_check_l && !new_check_l){
            territory_l=hunting(jumping_place_l);
            if(!(find(territory.begin(),territory.end(),territory_l)!=territory.end())){
                territory.push_back(territory_l);
                final_access=access;
                new_check_l=true;
            }
        }
        //cout<<territory_l<<endl;
        count++;
    }while ((territory.size()<3)&&(territory.size()!=2||jumping_place_l>=low_bound)&&(territory.size()!=2||jumping_place_r<=up_bound)&&(territory.size()!=1||jumping_place_r<=up_bound||jumping_place_l>=low_bound));
    int dir=0;
    float territory_distance;
    bool reverse_tr= false;
    if(best(territory_l,now_territory,territory_r)==territory_l){
        dir=-1;
        if((territory_l-now_territory)>0){
            reverse_tr=true;
            dir=1;
        }
    }
    else if(best(territory_l,now_territory,territory_r)==now_territory){
        int counter=0;
        bool jumping_speed_down_check= false;
        do{
            //cout<<now_territory<<"====="<<endl;
            //cout<<"--------------------"<<endl;
            counter--;
            float next_place_r,next_place_l;
            float next_territory_r,next_territory_l;
            next_place_r=now_territory+pow(2,counter)*this_time_radius;
            next_place_l=now_territory-pow(2,counter)*this_time_radius;
            if(next_place_r>up_bound)
                next_place_r=up_bound;
            if(next_place_l<low_bound)
                next_place_l=low_bound;
            epin_radius=pow(2,counter)*this_time_radius;
            next_territory_r=hunting(next_place_r);
            if(!(find(territory.begin(),territory.end(),next_territory_r)!=territory.end())){
                territory.push_back(next_territory_r);
                final_access=access;
            }
            epin_radius=-(pow(2,counter)*this_time_radius);
            next_territory_l=hunting(next_place_l);
            //cout<<next_territory_r<<endl;
            //cout<<next_territory_l<<endl;
            if(!(find(territory.begin(),territory.end(),next_territory_l)!=territory.end())){
                territory.push_back(next_territory_l);
                final_access=access;
            }
            if(next_territory_r==next_territory_l) {
                jumping_speed_down_check = true;
            }
            now_territory=min(territory);
        }while(!jumping_speed_down_check);
    }
    else if(best(territory_l,now_territory,territory_r)==territory_r){
        dir=1;
        if((territory_r-now_territory)<0){
            reverse_tr=true;
            dir=-1;
        }
    }
    if(dir){
        if(dir==1){
            if(reverse_tr==true){
                territory_distance=fabs(now_territory-territory_l);
                now_territory=territory_l;
            }
            else{
                    territory_distance = fabs(now_territory - territory_r);
                    now_territory=territory_r;
            }
        }
        if(dir==-1){
            if(reverse_tr==true){
                territory_distance=fabs(now_territory-territory_r);
                now_territory=territory_r;
            }
            else {
                territory_distance = fabs(now_territory - territory_l);
                now_territory = territory_l;
            }
        }
        float next_place;
        float next_territory;
        int counter=1;
        next_place=now_territory+pow(2,counter)*dir*territory_distance;
        if(next_place>up_bound)
            next_place=up_bound;
        if(next_place<low_bound)
            next_place=low_bound;
        epin_radius=pow(2,counter)*dir*territory_distance;
        next_territory=hunting(next_place);
        //cout<<next_territory<<endl;
        if(!(find(territory.begin(),territory.end(),next_territory)!=territory.end())){
            territory.push_back(next_territory);
            final_access=access;
        }
        do{
            if(function(now_territory)>function(next_territory)){
                counter++;
                now_territory=next_territory;
                next_place=now_territory+pow(2,counter)*dir*territory_distance;
                if(next_place>up_bound)
                    next_place=up_bound;
                if(next_place<low_bound)
                    next_place=low_bound;
                epin_radius=pow(2,counter)*dir*territory_distance;
                next_territory=hunting(next_place);
                //cout<<next_territory<<endl;
                if(!(find(territory.begin(),territory.end(),next_territory)!=territory.end())){
                    territory.push_back(next_territory);
                    final_access=access;
                }
            }

        }while(function(now_territory)>function(next_territory));
        bool jumping_speed_down_check= false;
        do{
            //cout<<now_territory<<"====="<<endl;
            //cout<<"--------------------"<<endl;
            counter--;
            float next_place_r,next_place_l;
            float next_territory_r,next_territory_l;
            next_place_r=now_territory+pow(2,counter)*territory_distance;
            next_place_l=now_territory-pow(2,counter)*territory_distance;
            if(next_place_r>up_bound)
                next_place_r=up_bound;
            if(next_place_l<low_bound)
                next_place_l=low_bound;
            epin_radius=pow(2,counter)*territory_distance;
            next_territory_r=hunting(next_place_r);
            if(!(find(territory.begin(),territory.end(),next_territory_r)!=territory.end())){
                territory.push_back(next_territory_r);
                final_access=access;
            }
            epin_radius=-(pow(2,counter)*territory_distance);
            next_territory_l=hunting(next_place_l);
            //cout<<next_territory_r<<endl;
            //cout<<next_territory_l<<endl;
            if(!(find(territory.begin(),territory.end(),next_territory_l)!=territory.end())){
                territory.push_back(next_territory_l);
                final_access=access;
            }
            if(next_territory_r==next_territory_l) {
                jumping_speed_down_check = true;
            }
            now_territory=min(territory);
        }while(!jumping_speed_down_check);
    }


}
// Let variable data Initialize in Parameter instead of in a constant function
// Ctrl + Alt + L would reformat your code just like gofmt, try it.


int main() {
    f << "Dimension : 1" << endl;
    f<<"Formula : 10+sum(Xi^2-10*cos(2*Pi*Xi))"<<endl;
    f<<"Range : -15 ~ 15"<<endl;
    f<<"Position :"<<endl;
    srand(114);
    Jaguar jaguar;
    jaguar.access=0;
    for(int i=0;i<30;i++) {
        jaguar.access=0;
        jaguar.hunting();
        jaguar.jumping();
        cout<<jaguar.final_access<<endl;
        //for(int i=0;i<jaguar.territory.size();i++){
            //cout<<setprecision(10)<<jaguar.territory.at(i)<<"   ";
        //}
        jaguar.territory.clear();
    }
    f.close();
    return 0;
}