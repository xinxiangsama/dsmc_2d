# Compiler and flags
# 定义编译器
CC = mpicxx

# 定义编译选项，这里可以根据需要添加优化选项等
# CCFLAGS = -O3 -march=native -funroll-loops -ffast-math -flto -Wunused-result -Wno-return-type -fprefetch-loop-arrays -I./src -I./include -I/home/xinxiangsama/lib/hdf5-1.14.6/include -DH5_BUILD_AS_DYNAMIC_LIB
CCFLAGS = -O3 -Wunused-result -Wno-return-type -I./src -I./include -I/home/xinxiangsama/lib/hdf5-1.14.6/include -DH5_BUILD_AS_DYNAMIC_LIB

# 定义链接选项
# LDFLAGS = -lboost_system -lboost_filesystem -L/home/xinxiangsama/lib/hdf5-1.14.6/lib -lhdf5_cpp -lhdf5
LDFLAGS = -L/home/xinxiangsama/lib/hdf5-1.14.6/lib -lhdf5_cpp -lhdf5

# 定义源文件
SOURCES = \
	./src/main.cpp \
	./src/Run.cpp \
	./src/particle/Particle.cpp \
	./src/cell/Cell.cpp \
	./src/meshes/Element.cpp  ./src/meshes/CartesianMesh.cpp \
	./src/parallel/Parallel.cpp ./src/parallel/CartesianParallel.cpp \
	./src/random/Random.cpp \
	./src/boundary/WallBoundary.cpp ./src/boundary/OutflowBoundary.cpp ./src/boundary/PeriodicBoundary.cpp\
	./src/io/Output.cpp \
	./src/object/Segment.cpp ./src/object/Circle.cpp \
	

# 定义目标文件
OBJECTS = $(patsubst %.cpp, build/%.o, $(SOURCES))

# 定义最终的可执行文件名
EXECUTABLE = runsim

# 默认目标
all: $(EXECUTABLE)

# 链接目标文件生成可执行文件
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CCFLAGS) -o $@ $^ $(LDFLAGS)

# 从.cpp文件编译出.o文件
build/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CC) $(CCFLAGS) -c $< -o $@

# 清理规则
clean:
	rm -f $(OBJECTS) $(EXECUTABLE)
	rm -rf ./res/*.h5

# 防止make自动生成*.o文件的规则
.PHONY: all clean