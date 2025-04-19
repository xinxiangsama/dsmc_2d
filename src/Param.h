#pragma once
#include <cmath>

//==================物理参数=======================
static constexpr double boltz = 1.3806e-23; // J/K 波尔兹曼常数
static constexpr double mass = 6.63e-26; // mass argon 分子质量
static constexpr double diam = 4.17e-10; // eff diam argon 分子直径
static constexpr double Volume_Particle = M_PI * diam * diam; //分子的体积
static constexpr unsigned int N_Particle =6.4e7; //总模拟分子数
static constexpr double tau = 5e-8; //时间步长（s） 需要小于分子平均碰撞频率
static constexpr double Vtl = 0.81;
static constexpr double VSS_coe = 1.0; 

//=================流场参数========================
static constexpr double L1 = 16.0e-2; //横向维度 (m)
static constexpr double L2 = 16.0e-2; //纵向维度 (m)
static constexpr double L3 = 4.0e-4; 
static constexpr double Rho= 1.0e-3 ; // 流场密度 kg/m^3 
static constexpr double T = 273.0; // temperature (K)

//=================圆柱参数========================
static constexpr double Radius = 1.85e-2; //圆柱半径
static constexpr double Center_x = 8.0e-2; //圆柱中心横坐标
static constexpr double Center_y = 8.0e-2; //圆柱中心纵坐标
static constexpr double V_jet = 1300.0;//来流速度
//===============正交结构网格==================
static constexpr unsigned int N1 = 400; //横向网格数
static constexpr unsigned int N2 = 400; // 纵向网格数
static constexpr unsigned int N3 = 1; // 纵向网格数

//=============一些constexper=========
static constexpr double Volume = L1 * L2 * L3; //流场体积
static constexpr size_t Fn = (Volume * Rho / mass) / N_Particle; //每个仿真分子代表的物理分子数


//==============写出文件相关参数============
static constexpr int skip_every = 10; // skip every _ steps when writing out
static constexpr int n_zero = 5; // 10 digits for indexing output CSV series
static const int average_step = 10; // average over _ steps when writing out

