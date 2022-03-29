#Rm des précédantes itération
rm mpi1 omp rma rma_omp master slave *.txt

#Compilation des scripts
mpic++ -o mpi1 -I ../common/ ../vanilla.cpp ../common/fonctions.cpp 
mpic++ -o omp -fopenmp -I ../common/ ../omp.cpp ../common/fonctions.cpp  
mpic++ -o rma -I ../common/ ../rma.cpp ../common/fonctions.cpp 
mpic++ -o rma_omp -fopenmp -I ../common/ ../rma_omp.cpp ../common/fonctions.cpp
mpic++ -o master -I ../common/ ../master_slave/master.cpp ../common/fonctions.cpp 
mpic++ -o slave -I ../common/ ../master_slave/slave.cpp ../common/fonctions.cpp 

N=256
M=16

echo "Running with N=$N and M=$M\n"
for i in {1..10}
do
    echo "############ Iteration $i ############\n"
    
    echo "	mpi"
    mpirun -hostfile hostfile mpi1 $N $M 0 res.txt >> "mpi1-$N-$M.txt" 
    
    echo "	omp"
    mpirun -hostfile hostfile omp $N $M 0 res.txt >> "omp-$N-$M.txt" 
    
    echo "	rma"
    mpirun -hostfile hostfile rma $N $M 0 res.txt >> "rma-$N-$M.txt" 
    
    echo "	rma_omp"
    mpirun -hostfile hostfile rma_omp $N $M 0 res.txt >> "rma-omp-$N-$M.txt" 
    
    echo "	master"
    mpirun -np 1 -hostfile hostfile master $N $M 2 0 8 res.txt slave >> "ms-$N-$M.txt" 
done

N=1024
M=16

echo "Running with N=$N and M=$M\n"
for i in {1..10}
do
    echo "############ Iteration $i ############\n"

    echo "	mpi"
    mpirun -hostfile hostfile mpi1 $N $M 0 res.txt >> "mpi1-$N-$M.txt" 
    
    echo "	omp"
    mpirun -hostfile hostfile omp $N $M 0 res.txt >> "omp-$N-$M.txt" 
    
    echo "	rma"
    mpirun -hostfile hostfile rma $N $M 0 res.txt >> "rma-$N-$M.txt" 
    
    echo "	rma_omp"
    mpirun -hostfile hostfile rma_omp $N $M 0 res.txt >> "rma-omp-$N-$M.txt" 
    
    echo "	master"
    mpirun -np 1 -hostfile hostfile master $N $M 2 0 8 res.txt slave >> "ms-$N-$M.txt" 
done

N=4096
M=16

echo "Running with N=$N and M=$M\n"
for i in {1..10}
do
    echo "############ Iteration $i ############\n"

    echo "	mpi"
    mpirun -hostfile hostfile mpi1 $N $M 0 res.txt >> "mpi1-$N-$M.txt" 
    
    echo "	omp"
    mpirun -hostfile hostfile omp $N $M 0 res.txt >> "omp-$N-$M.txt" 
    
    echo "	rma"
    mpirun -hostfile hostfile rma $N $M 0 res.txt >> "rma-$N-$M.txt" 
    
    echo "	rma_omp"
    mpirun -hostfile hostfile rma_omp $N $M 0 res.txt >> "rma-omp-$N-$M.txt" 
    
    echo "	master"
    mpirun -np 1 -hostfile hostfile master $N $M 2 0 8 res.txt slave >> "ms-$N-$M.txt" 
done

N=256
M=256

echo "Running with N=$N and M=$M\n"
for i in {1..10}
do
    echo "############ Iteration $i ############\n"

    echo "	mpi"
    mpirun -hostfile hostfile mpi1 $N $M 0 res.txt >> "mpi1-$N-$M.txt" 
    
    echo "	omp"
    mpirun -hostfile hostfile omp $N $M 0 res.txt >> "omp-$N-$M.txt" 
    
    echo "	rma"
    mpirun -hostfile hostfile rma $N $M 0 res.txt >> "rma-$N-$M.txt" 
    
    echo "	rma_omp"
    mpirun -hostfile hostfile rma_omp $N $M 0 res.txt >> "rma-omp-$N-$M.txt" 
    
    echo "	master"
    mpirun -np 1 -hostfile hostfile master $N $M 2 0 8 res.txt slave >> "ms-$N-$M.txt" 
done

N=1024
M=256

echo "Running with N=$N and M=$M\n"
for i in {1..10}
do
    echo "############ Iteration $i ############\n"

    echo "	mpi"
    mpirun -hostfile hostfile mpi1 $N $M 0 res.txt >> "mpi1-$N-$M.txt" 
    
    echo "	omp"
    mpirun -hostfile hostfile omp $N $M 0 res.txt >> "omp-$N-$M.txt" 
    
    echo "	rma"
    mpirun -hostfile hostfile rma $N $M 0 res.txt >> "rma-$N-$M.txt" 
    
    echo "	rma_omp"
    mpirun -hostfile hostfile rma_omp $N $M 0 res.txt >> "rma-omp-$N-$M.txt" 
    
    echo "	master"
    mpirun -np 1 -hostfile hostfile master $N $M 2 0 8 res.txt slave >> "ms-$N-$M.txt" 
done

N=4096
M=256

echo "Running with N=$N and M=$M\n"
for i in {1..10}
do
    echo "############ Iteration $i ############\n"

    echo "	mpi"
    mpirun -hostfile hostfile mpi1 $N $M 0 res.txt >> "mpi1-$N-$M.txt" 
    
    echo "	omp"
    mpirun -hostfile hostfile omp $N $M 0 res.txt >> "omp-$N-$M.txt" 
    
    echo "	rma"
    mpirun -hostfile hostfile rma $N $M 0 res.txt >> "rma-$N-$M.txt" 
    
    echo "	rma_omp"
    mpirun -hostfile hostfile rma_omp $N $M 0 res.txt >> "rma-omp-$N-$M.txt" 
    
    echo "	master"
    mpirun -np 1 -hostfile hostfile master $N $M 2 0 8 res.txt slave >> "ms-$N-$M.txt" 
done
