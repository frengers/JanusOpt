import numpy as np
import time
from scoop import futures
import scoop
from deap import base, creator, tools
from deap import cma

import shutil, glob

import objective_function as fx

def test_func(x):
    cmd_list = []
    cmd_list.append('--tauc='+str(x[0]))
    cmd_list.append('--taucwepp='+str(x[1]))
    val = fx.run(cmd_list)
    return val['result'],
    

creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
creator.create("Individual", list, fitness=creator.FitnessMin)

toolbox = base.Toolbox()
toolbox.register("evaluate", test_func)

toolbox.register("map", futures.map)
strategy = cma.Strategy(centroid=[100, 20], sigma=25, lambda_=24)
toolbox.register("generate", strategy.generate, creator.Individual)
toolbox.register("update", strategy.update)

def main():
    np.random.seed(128)
    hof = tools.HallOfFame(1)
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stop = False
    for gen in range(10):
        # Generate a new population
        population = toolbox.generate()
        # Evaluate the individuals
        fitnesses = toolbox.map(toolbox.evaluate, population)
        for ind, fit in zip(population, fitnesses):
            print fit
            ind.fitness.values = fit
        # Update the strategy with the evaluated individuals
        toolbox.update(population)
        hof.update(population)
        stats.update(population)
        
        # Delete directories
        trials = glob.glob('trial-*')
        for t in trials:
            shutil.rmtree(t)
            
        print hof    
        
    return hof


if __name__ == '__main__':
    
    x = main()
    print time.ctime(), "finished", x


    

    
