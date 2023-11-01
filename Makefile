
CXX=g++
CFLAGS = -std=gnu++17 -Wall -fopenmp -O2 -w
INC=-I include/

index: clean_index
	$(CXX) $(INC) $(CFLAGS) src/MurmurHash3.cpp src/Rambo_constructionIDL.cpp src/bitArray.cpp src/utils.cpp src/MyBloom.cpp src/IDLBloomFilter.cpp \
							src/main.cpp \
	-o index $(LFLAGS)

query: clean_query
	$(CXX) $(INC) $(CFLAGS) src/MurmurHash3.cpp src/Rambo_constructionIDL.cpp src/bitArray.cpp src/utils.cpp src/MyBloom.cpp src/IDLBloomFilter.cpp \
							src/query.cpp \
	-o query $(LFLAGS)
	
clean_index:
	rm -f index
clean_query:
	rm -f query

.PHONY: clean all

debug: CXXFLAGS += -DDEBUG -g
debug: all

release: CXXFLAGS += -O2
release: all
