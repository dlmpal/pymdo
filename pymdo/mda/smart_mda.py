from typing import List
from typing import Dict

import networkx as nx

from pymdo.core.discipline import Discipline
from .mda import MDA
from .factory import mda_factory


class SmartMDA(MDA):
    """
    
    This MDA sub-class creates an execution sequence for the disciplines,
    according to their inter-depedencies.

    This can reduce the execution time, by excluding weakly-coupled disciplines from
    MDA loops.

    By default, Gauss-Seidel MDAs are created to resolve couplings that might emerge. 
    The user can also specify which MDA algorithm to use.
    
    """

    def __init__(self,
                 disciplines: List[Discipline],
                 name: str = "SmartMDA"):

        super().__init__(disciplines,
                         name)

        self.groups: List[Dict[str, MDA]] = self._CreateGroups()

    def _CreateGroups(self) -> List[Dict[str, MDA]]:
        """

        Create discipline groups.

        Groups are split into levels. The groups in each level can execute in parallel.

        Each group is either a single discipline, or an MDA of coupled disciplines.

        """

        self.groups = []

        graph = nx.DiGraph()

        for disc_i in self.disciplines:

            for var_i in disc_i.outputVars:

                for disc_j in self.disciplines:

                    if var_i in disc_j.inputVars:

                        graph.add_edge(disc_i,
                                       disc_j)

        groupList: List[List[Discipline]] = []

        for group in sorted(nx.strongly_connected_components(graph)):

            groupList.append(list(group))

        graphCondensed = nx.condensation(
            graph, nx.strongly_connected_components(graph))

        execSeq = []

        while True:

            if len(graphCondensed.nodes) == 0:

                break

            currentLevel = []

            for groupIdx in graphCondensed:

                if graphCondensed.out_degree[groupIdx] == 0:

                    currentLevel.append(groupIdx)

            execSeq.append(currentLevel)

            for groupIdx in currentLevel:

                graphCondensed.remove_node(groupIdx)

        for level in execSeq[::-1]:

            curLevelGroups = {}

            for groupIdx in level:

                groupDisciplines = groupList[groupIdx]

                groupName = f"Group_{groupIdx}"

                
                if len(groupDisciplines) > 1:
                    # Group contains strongly coupled disciplines,
                    # therefore an mda is used
                    curLevelGroups[groupName] = mda_factory(groupDisciplines)
                
                else:
                    # Group is just a single discipline
                    # No mda needed 
                    curLevelGroups[groupName] = groupDisciplines[0]

            self.groups.append(curLevelGroups)

        return self.groups

    def _eval(self) -> None:

        for lvl in self.groups:

            for group in lvl.values():

                self.values.update(group.eval(self.values))
