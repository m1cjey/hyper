// stdafx.h : 標準のシステム インクルード ファイルのインクルード ファイル、または
// 参照回数が多く、かつあまり変更されない、プロジェクト専用のインクルード ファイル
// を記述します。
//




#pragma once

#include "targetver.h"

#include <stdio.h>
#include <tchar.h>
#include <time.h>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <Windows.h>
#include <fstream>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <algorithm>
#include <string>
#include <sstream>
#include <omp.h>


#include "tetgen.h"
#include "define.h"
#include "CONFIG.h"
#include "PART.h"
#include "FEM3D.h"
#include "tetgen_config.h"
#include "Hyperelastic.h"
#include "tetfunc.h"
#include "elastic_config.h"
#include "function.h"
#include "elastic_calculation.h"
#include "math_library.h"

// TODO: プログラムに必要な追加ヘッダーをここで参照してください。

#define _CRTDBG_MAP_ALLOC//メモリリーク検出用
#include <stdlib.h>
#include <crtdbg.h>//メモリリーク検出用
#define new ::new(_NORMAL_BLOCK, __FILE__, __LINE__)//メモリリーク検出用





