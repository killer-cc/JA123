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
int dim=1;
int dim_count;
// Use a function pointer instead of constant function for extendability
// Example: std::function<double(float)> FitnessFunction;


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
    vector<vector<float>> territory;
    const float up_bound = 15.0;
    const float low_bound = -15.0;
    float punish = 2147483647.0;
    vector<float> place;
    double fitness;
    vector<float> past_best_place;
    double past_best_fitness;
    float step;
    bool right_or_left;
    bool hunt_check ;
    int acc_step_count;
    float access;
    float jumping_radius;
    vector<float> now_territory;
    float epin_radius;
    int final_access;


    Jaguar();

    void init();

    void init(float);

    void look_around();

    void hunting_jump();

    vector<float> speed_down_around();

    void hunting();

    vector<float> hunting(float);

    void jumping();

    double function(float);

    float best(float,float,float);

    vector<float> min(vector<vector<float>>);
};
Jaguar::Jaguar() {
    for(int i=0;i<dim;i++){
        place.push_back(0);
    }
    fitness = 0;
}
double Jaguar::function(float a) {
    //return 20.0+double(exp(1))-20.0*double(exp(-0.2*double(sqrt(double(pow(a,2.0))))))-double(exp(double(cos(2.0*double(acos(-1.0))*a))));
    /*if(true){
        double b=0,c=0,d=0;
        for(int i=0;i<dim;i++){
            if (i!= dim_count){
                b=b+double(pow(place[i],2.0));
            } else {
                b=b+double(pow(a,2.0));
            }
        }
        b=b/double(dim);
        b=double(sqrt(b));
        for(int i=0;i<dim;i++){
            if (i!= dim_count){
                c=c+double(cos(2.0*double(acos(-1.0))*place[i]));
            } else {
                c=c+double(cos(2.0*double(acos(-1.0))*a));
            }
        }
        c=c/double(dim);
        d=-20.0*double(exp(-0.2*b))-double(exp(c))+double(exp(1))+20.0;
        //be aware the order
        return d;
    }*/
    //return 10.0 + (pow(a, 2.0) - (10.0 * double(cos(2.0 * double(acos(-1.0)) * a))));
    /*if(true){
        double b=0;
        for(int i=0;i<dim;i++){
            if (i!= dim_count){
                b=b+(pow(place[i], 2.0) - (10.0 * (cos(2.0 * (acos(-1.0)) * place[i]))));
            } else {
                b=b+(pow(a, 2.0) - (10.0 * (cos(2.0 * (acos(-1.0)) * a))));
            }
        }
        b=b+10.0*dim;
        //be aware the order
        return b;
    }*/
    if (true){
        double b=0;
        for (int i=0;i<dim;i++){
            if(i!=dim_count){
                b=b+fabs(place[i]);
            } else{
                b=b+fabs(a);
            }
        }
    }
    /*if (true){
        double b=0;
        for (int i=0;i<dim;i++){
            if(i!=dim_count){
                b=b+(place[i]*place[i]);
            } else{
                b=b+(a*a);
            }
        }
    }*/
}
float Jaguar::best(float a, float b, float c){
    float best;
    best=float(function(a))<float(function(b))?a:b;
    best=float(function(best))<float(function(c))?best:c;
    //float to compare
    return best;
}
vector<float> Jaguar::min(vector<vector<float>> a){
    vector<float> min=a.at(0);
    for(int i=0;i<a.size();i++){
        if(float(function(a.at(i)[dim_count]))<float(function(min[dim_count])))
            //float to compare
            min=a.at(i);
    }
    return min;
}
void Jaguar::init() {
    right_or_left = false;
    hunt_check = false;
    acc_step_count = 1;
    jumping_radius=0;
    if(dim_count == 0){
        for (int i = 0; i < dim; i++) {
            place[i] = rand() / float(RAND_MAX) * (up_bound - low_bound) + low_bound;
        }
    }
    /*for(int i=0;i<dim;i++){
        cout << setprecision(20) << place[i]<<" ";
    }
     cout<< endl;*/
    jumping_radius=place[dim_count];
    past_best_place = place;
    fitness = function(place[dim_count]);
    //f<<setprecision(20)<<"*"<<epin_count<<" "<<fitness<<": "<<place<<",10,0,12"<<endl;
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
    place[dim_count] = a;
    /*for(int i=0;i<dim;i++){
        cout << setprecision(20) << place[i]<<" ";
    }
    cout<< endl;*/
    past_best_place = place;
    jumping_radius=place[dim_count];
    fitness = function(place[dim_count]);
    //f<<setprecision(20)<<"*"<<epin_count<<" "<<fitness<<": "<<place<<",10,"<<epin_radius<<",12"<<endl;
    epin_count++;
    past_best_fitness = fitness;
    access++;
    step = powf(2.0, int(log2(up_bound)) - 11);
}
void Jaguar::look_around() {
    float place_l, place_r;
    double fitness_l, fitness_r;
    place_l = place[dim_count] - step * acc_step_count;
    place_r = place[dim_count] + step * acc_step_count;
    if (place_r < low_bound || place_r > up_bound)
        fitness_r = punish;
    else{
        fitness_r = function(place_r);
        access++;
    }
    /*for(int i=0;i<dim;i++){
        if(i==dim_count){
            cout<< setprecision(20)<<place_r<<" ";
        }else{
            cout<< setprecision(20)<<place[i]<<" ";
        }
    }
    cout<< setprecision(20)<<"," << step * acc_step_count<< "," << access << ","<< fitness_r << endl;*/
    f<<setprecision(20)<<"*"<<epin_count<<" "<<fitness_r<<": "<<place_r<<",10,0,12"<<endl;
    epin_count++;
    if (place_l < low_bound || place_l > up_bound)
        fitness_l = punish;
    else{
        fitness_l = function(place_l);
        access++;
    }
    /*for(int i=0;i<dim;i++){
        if(i==dim_count){
            cout<< setprecision(20)<<place_l<<" ";
        }else{
            cout<< setprecision(20)<<place[i]<<" ";
        }
    }
    cout<< setprecision(20)<<"," << step * acc_step_count<< "," << access << ","<< fitness_l << endl;*/
    f<<setprecision(20)<<"*"<<epin_count<<" "<<fitness_l<<": "<<place_l<<",10,0,12"<<endl;
    epin_count++;
    place[dim_count] = fitness < fitness_r ? place[dim_count] : place_r;
    fitness = fitness < fitness_r ? fitness : fitness_r;
    place[dim_count]=fitness<fitness_l?place[dim_count]:place_l;
    fitness = fitness<fitness_l?fitness:fitness_l;
    if (past_best_fitness>fitness  ) {
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
        vr_place = place[dim_count] + step * acc_step_count*dir;
        if (vr_place > low_bound && vr_place < up_bound) {
            vr_fitness = function(vr_place);
            access++;
        } else {
            vr_fitness = punish;
        }
        /*for(int i=0;i<dim;i++){
            if(i==dim_count){
                cout<< setprecision(20)<<vr_place<<" ";
            }else{
                cout<< setprecision(20)<<place[i]<<" ";
            }
        }
        cout<< setprecision(20)<<"," << step * acc_step_count<< "," << access << ","<< vr_fitness << endl;*/
        f<<setprecision(20)<<"*"<<epin_count<<" "<<vr_fitness<<": "<<vr_place<<",10,0,12"<<endl;
        epin_count++;
        if (vr_fitness < fitness) {
            place[dim_count] = vr_place;
            fitness = vr_fitness;
            past_best_place = place;
            past_best_fitness = fitness;
        } else {
            acc_step_count = acc_step_count / 2;
            vr_place = place[dim_count] + step * acc_step_count*dir;
            if (vr_place > low_bound && vr_place < up_bound) {
                vr_fitness = function(vr_place);
                access++;
            } else {
                vr_fitness = punish;
            }
            place[dim_count] = vr_place;
            fitness = vr_fitness;
            if (past_best_fitness > fitness) {
                past_best_place = place;
                past_best_fitness = fitness;
            }
            //cout << setprecision(20) << place[dim_count] << "," << step * acc_step_count << "," << access<< "," << fitness << endl;
            /*for(int i=0;i<dim;i++){
                cout<< setprecision(20)<<place[i]<<" ";
            }
            cout<< setprecision(20)<<"," << step * acc_step_count<< "," << access << ","<< fitness << endl;*/
            //f<<setprecision(20)<<"*"<<epin_count<<" "<<fitness<<": "<<place<<",10,0,12"<<endl;
            epin_count++;
            stop_flag = true;
        }
    }
}
vector<float> Jaguar::speed_down_around() {
    while (acc_step_count > 1) {
        acc_step_count = acc_step_count / 2;
        place = past_best_place;
        fitness = past_best_fitness;
        look_around();
    }
    double a = past_best_fitness;
    while (past_best_fitness == a) {
        place=past_best_place;
        fitness=past_best_fitness;
        step = step / 2;
        if ((place[dim_count]+step)==place[dim_count]||(place[dim_count]-step)==place[dim_count]) {
            //one side to lose precision and stop
            jumping_radius=fabs(jumping_radius-place[dim_count]);
            hunt_check=true;
            return place;
        }
        else {
            look_around();
        }
    }
    return place;
}
void Jaguar::hunting() {
    init();
    float a=place[dim_count];
    look_around();
    vector<float> final_place;
    while(a==place[dim_count]){
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
vector<float> Jaguar::hunting(float a) {
    hunt_check= false;
    vector<float> final_place;
    init(a);
    float b=place[dim_count];
    look_around();
    while(b==place[dim_count]){
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
    now_territory=territory.at(territory.size()-1);
    int now_territory_place=territory.size()-1;
    final_access=access;
    vector<float> territory_r,territory_l;
    int count=0;
    float jumping_place_r,jumping_place_l;
    float this_time_radius=jumping_radius;
    float this_time_radius_mid=this_time_radius;
    //radius for mid is better
    bool outside_check_r= false,outside_check_l= false;
    bool new_check_r=false,new_check_l=false;
    bool left_find_tr=false,right_find_tr=false;
    float sw_this_time_radius_r=this_time_radius*2.0f,sw_this_time_radius_l=this_time_radius*2.0f;
    do{

        if(count!=0){
            float a;
            //cout<<fabs(jumping_place_r-now_territory)<<endl;
            //cout<<fabs(jumping_place_l-now_territory)<<endl;
            if(!outside_check_r && !new_check_r&&!outside_check_l && !new_check_l&&fabs(jumping_place_r-now_territory[dim_count])>=fabs(jumping_place_l-now_territory[dim_count])){
                a=fabs(jumping_place_r-now_territory[dim_count]);
            }
            else{
                a=fabs(jumping_place_l-now_territory[dim_count]);
            }
            if(outside_check_r || new_check_r)
            {
                a=fabs(jumping_place_l-now_territory[dim_count]);
            }
            if(outside_check_l || new_check_l){
                a=fabs(jumping_place_r-now_territory[dim_count]);
            }
            this_time_radius_mid=a;
            sw_this_time_radius_r=a*2.0f;
            sw_this_time_radius_l=a*2.0f;
        }
        epin_radius=sw_this_time_radius_r;
        jumping_place_r=territory.at(now_territory_place)[dim_count]+sw_this_time_radius_r;
        if(!outside_check_r && jumping_place_r>=up_bound){
            territory_r=hunting(up_bound);
            if(!(find(territory.begin(),territory.end(),territory_r)!=territory.end())){
                territory.push_back(territory_r);
                final_access=access;
                new_check_r=true;
            }
            outside_check_r= true;
            right_find_tr=true;
        }
        if (!outside_check_r && !new_check_r){
            territory_r=hunting(jumping_place_r);
            if(!(find(territory.begin(),territory.end(),territory_r)!=territory.end())){
                territory.push_back(territory_r);
                final_access=access;
                new_check_r=true;
                right_find_tr=true;
            }
            if(territory_r!=now_territory){
                right_find_tr=true;
            }
        }
        //cout<<territory_r<<endl;
        epin_radius=sw_this_time_radius_l;
        jumping_place_l=territory.at(now_territory_place)[dim_count]-sw_this_time_radius_l;
        if(!outside_check_l && jumping_place_l<=low_bound){
            territory_l=hunting(low_bound);
            if(!(find(territory.begin(),territory.end(),territory_l)!=territory.end())){
                territory.push_back(territory_l);
                final_access=access;
                new_check_l=true;
            }
            outside_check_l= true;
            left_find_tr=true;
        }
        if (!outside_check_l && !new_check_l){
            territory_l=hunting(jumping_place_l);
            if(!(find(territory.begin(),territory.end(),territory_l)!=territory.end())){
                territory.push_back(territory_l);
                final_access=access;
                new_check_l=true;
                left_find_tr=true;
            }
            if(territory_l!=now_territory){
                left_find_tr=true;
            }
        }
        //cout<<territory_l<<endl;
        count++;
    }while (!left_find_tr || !right_find_tr);
    int dir=0;
    float territory_distance;
    bool reverse_tr= false;
    bool same_tr=false;
    //checker for three point is same and go into mid is better function
    if (territory_l[dim_count]==now_territory[dim_count] && now_territory[dim_count]==territory_r[dim_count]){
        same_tr=true;
    }
    if(best(territory_l[dim_count],now_territory[dim_count],territory_r[dim_count])==territory_l[dim_count] && !same_tr){
        dir=-1;
        if((territory_l[dim_count]-now_territory[dim_count])>0){
            reverse_tr=true;
            dir=1;
        }
    }
    else if(best(territory_l[dim_count],now_territory[dim_count],territory_r[dim_count])==now_territory[dim_count]){
        int counter=0;
        bool jumping_speed_down_check= false;
        do{
            //cout<<now_territory<<"====="<<endl;
            //cout<<"--------------------"<<endl;
            counter--;
            float next_place_r,next_place_l;
            vector<float> next_territory_r,next_territory_l;
            next_place_r=now_territory[dim_count]+pow(2,counter)*this_time_radius_mid;
            next_place_l=now_territory[dim_count]-pow(2,counter)*this_time_radius_mid;
            if(next_place_r>up_bound)
                next_place_r=up_bound;
            if(next_place_l<low_bound)
                next_place_l=low_bound;
            epin_radius=pow(2,counter)*this_time_radius_mid;
            next_territory_r=hunting(next_place_r);
            if(!(find(territory.begin(),territory.end(),next_territory_r)!=territory.end())){
                territory.push_back(next_territory_r);
                final_access=access;
            }
            epin_radius=-(pow(2,counter)*this_time_radius_mid);
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
            int cc=counter-1;
            if(((now_territory[dim_count]+float(pow(2,cc)*this_time_radius_mid))==now_territory[dim_count]) || ((now_territory[dim_count]-float(pow(2,cc)*this_time_radius_mid))==now_territory[dim_count])){
                jumping_speed_down_check=true;
                place=now_territory;
            }
        }while(!jumping_speed_down_check);
    }
    else if(best(territory_l[dim_count],now_territory[dim_count],territory_r[dim_count])==territory_r[dim_count]&& !same_tr){
        dir=1;
        if((territory_r[dim_count]-now_territory[dim_count])<0){
            reverse_tr=true;
            dir=-1;
        }
    }
    if(dir){
        if(dir==1){
            if(reverse_tr==true){
                territory_distance=fabs(now_territory[dim_count]-territory_l[dim_count]);
                now_territory=territory_l;
            }
            else{
                    territory_distance = fabs(now_territory[dim_count] - territory_r[dim_count]);
                    now_territory=territory_r;
            }
        }
        if(dir==-1){
            if(reverse_tr==true){
                territory_distance=fabs(now_territory[dim_count]-territory_r[dim_count]);
                now_territory=territory_r;
            }
            else {
                territory_distance = fabs(now_territory[dim_count] - territory_l[dim_count]);
                now_territory = territory_l;
            }
        }
        float next_place;
        vector<float> next_territory;
        int counter=1;
        next_place=now_territory[dim_count]+pow(2,counter)*dir*territory_distance;
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
            if(float(function(now_territory[dim_count]))>float(function(next_territory[dim_count]))){
                counter++;
                now_territory=next_territory;
                next_place=now_territory[dim_count]+pow(2,counter)*dir*territory_distance;
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

        }while(float(function(now_territory[dim_count]))>float(function(next_territory[dim_count])));
        bool jumping_speed_down_check= false;
        do{
            //cout<<now_territory<<"====="<<endl;
            //cout<<"--------------------"<<endl;
            counter--;
            float next_place_r,next_place_l;
            vector<float> next_territory_r,next_territory_l;
            next_place_r=now_territory[dim_count]+pow(2,counter)*territory_distance;
            next_place_l=now_territory[dim_count]-pow(2,counter)*territory_distance;
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
            if(next_territory_r[dim_count]==next_territory_l[dim_count]) {
                jumping_speed_down_check = true;
            }
            now_territory=min(territory);
            //checker for jumper is lose precision
            int cc=counter-1;
            if(((now_territory[dim_count]+float(pow(2,cc)*territory_distance))==now_territory[dim_count]) || ((now_territory[dim_count]-float(pow(2,cc)*territory_distance))==now_territory[dim_count])){
                jumping_speed_down_check=true;
                place=now_territory;
            }
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
    dim_count=0;
    for(int kappa=0;kappa<30;kappa++) {
        jaguar.access = 0;
        for (int i = 0; i < dim; i++) {
            jaguar.hunting();
            jaguar.jumping();
            //cout<<jaguar.final_access<<endl;
            //for(int i=0;i<jaguar.territory.size();i++){
            //cout<<setprecision(10)<<jaguar.territory.at(i)<<"   ";
            //}
            dim_count++;
        }
        dim_count=0;
        /*for (int i=0;i<jaguar.territory.size();i++){
            for(int j=0;j<dim;j++){
                cout<<setprecision(20)<<jaguar.territory[i][j]<<"  ";
            }
            cout<<endl;
        }*/
        cout<<jaguar.territory.size()<<endl;
        cout<<jaguar.final_access<<endl;
        jaguar.territory.clear();
        //cout<<jaguar.access<<endl;
    }
    f.close();
    return 0;
}