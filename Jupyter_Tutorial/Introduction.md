# Robustness and Plasticity in Regulatory Networks

* [Introduction](#Introduction)
    * [Robustness in biological systems](#Robustness-in-biological-systems)
    * [Boolean regulatory networks](#Boolean-Regulatory-Networks)
    * [Biological system: Th17/iTreg network](#Biological-system:-Th17/iTreg-network)
* [BoolNet](./RPRN-BoolNet.ipynb)
    * BoolNet and BoolNetPerturb
    * Construction
    * Attractors
    * Labels
* [Functions](./RPRN-Functions.ipynb)
    * Over-expression and KnockOuts
    * Fixed environments
    * Truth tables
* [Updating](./RPRN-Updating.ipynb)
    * Synchronous vs asynchronous
    * Transition table
* [States and trajectories](./RPRN-States-Trajectories.ipynb)
    * Transient perturbations
    * Stochastic noise
* [Appendix](./RPRN-Appendix.ipynb)
    * Simulation and model checking
    * Importing and exporting with SBML-qual



# Introduction

[//]: # (## About)

This tutorial explains how to use [BoolNet](https://cran.r-project.org/web/packages/BoolNet/index.html) and [BoolNetPerturb](https://github.com/mar-esther23/boolnet-perturb) to study the robustness of Boolean regulatory networks and the biological implications of this networks.

This tutorial supposes that the reader is familiar with the basic concepts of:
* __Boolean regulatory networks__. There are a lot of basic introductions to the topic like [Kaplan & Glass, 1995](https://books.google.com.mx/books?id=kfmThocv1qsC&pg=PA55#v=onepage&q&f=false), [Azpeitia 2011](http://journal.frontiersin.org/article/10.3389/fpls.2011.00092/abstract) and [Albert 2014](http://onlinelibrary.wiley.com/doi/10.1002/wsbm.1273/full). [CoLoMoTo](http://www.colomoto.org/) also published a more advanced review [CoLoMoTo 2015](http://www.biorxiv.org/content/biorxiv/early/2014/10/19/010504.full.pdf) full of useful references. 
* __R programming language__. A good starting point is the [R programming course](https://www.coursera.org/course/rprog) at Coursera.
* __Molecular biology__. Contact your local biologist.
* __Robustness__. The book [Wagner 2005](https://books.google.com.mx/books?id=pRFYAQAAQBAJ&printsec=frontcover#v=onepage&q&f=false) was a great inspiration for this work, the introduction is a good place to start.



## Robustness and plasticity in biological systems

Organisms develop in a changing world where they are subject to both intrinsic and extrinsic perturbations. They need to be both resilient and capable of altering their phenotype depending of the circumstances. This two behaviors coexist in all organisms, which suggests that there are common mechanism that underlies both robustness and plasticity. 

Robustness is the capacity of an organism of maintaining its biological function in response to perturbations, while plasticity is an organism's capacity of changing from one state or function to an other in response to perturbations. However, for studying robustness and plasticity it is necessary to determine _what_ function of the system is robust and to _which_ kind of perturbations.

Boolean Regulatory Networks (RN) are a useful tool for studying the robustness and plasticity of biological systems in response to different kinds of perturbations. 
RN integrate the molecular regulation to predict cellular level phenomena using a mathematical formalism. The model recovers the attractors of the system, that correspond to the biological cell types. 
The robustness and plasticity of this attractors in response to perturbations of the model can be used to study the robustness and plasticity of the biological system. 

We can say that an _attractor_ is robust to a perturbation if the RN returns to the same attractor, or plastic to a perturbation if the RN changes to a different attractor. 
The robustness or plasticity of the _system_ is the result of the individual robustness and plasticity of all its attractors. In this way, robustness is a characteristic of the system that emerges from the dynamic properties of its elements in response to perturbations.

The robustness of the attractors and of the system will be differ according to the perturbation. RN let us study multiple kinds of perturbations depending in which part of the RN we alter.



## Boolean Regulatory Networks

* The different levels of RN perturbations
* Attractors as initial conditions
* Modeling perturbations in GRN and biological implications


Boolean regulatory Networks (RN) are deterministic dynamic systems. RN have been used to study the differentiation, robustness and plasticity of developmental processes in different organisms. 
GRN consist of nodes -that represent genes, proteins, biological processes- and edges -that represent the regulatory interactions. Using this information it is possible to construct functions that describe the state of the nodes, this means, whether the genes or proteins of the system are active or inactive. 
The functions of the network are iterated to obtain the stable states or the system or attractors. This attractors correspond to the biological cell types.

<figure>
    <a href="images/DiagramGRN.png">
       <img src="images/DiagramGRN.png" alt="Boolean regulatory network.">
    </a>
<figcaption>
    Boolean regulatory networks.
</figcaption>
</figure>

However, as we have already discussed, biological systems are subjected to noise and perturbations. This perturbations can affect each of the levels of the network and correspond to different biological phenomena.



## Biological system: Th17/iTreg network

In this work we will study the Th17/iTreg regulatory network. CD4 + T cells are fundamental for the adaptive immune response. They integrate the signals of the environment and differentiate into different cell types (Th1, Th2, Th17, iTreg, etc), which activate different parts of the immune system. In particular, Th17 cells have been associated with the inflammatory response and iTreg cells with the regulation of the inflammatory response.

CD4+ T cells begin as naïve Th0 cells, which do not express a transcription factor. These cells are activated by antigen presentation and differentiate into different cell types depending in the cytokines in the environment. In the presence of IL-6 or IL-21 and TGF$\beta$ Th0 cells differentiate into Th17 cells and express ROR$\gamma$t, IL-21 and IL-17. In the presence of IL-2 and TGF$\beta$ Th0 cells differentiate into iTreg cells and express Foxp3 and TGF$\beta$. Interleukin-10 (IL-10) is a regulatory cytokine that is expressed in multiple cell types, it is activated by various cytokines including IL-6, IL-21 and TGF$\beta$. These cytokines and transcription factors regulate each other and their relationships can be visualized a as graph.

In this tutorial we will study the robustness and plasticity of the Th17/iTreg regulatory network in response to various kinds of perturbations.  


<figure>
    <a href="images/minTh17iTreg.png">
       <img src="images/minTh17iTreg.png" alt="Th17/iTreg regulatory network">
    </a>
<figcaption>
    Th17/iTreg regulatory network
</figcaption>
</figure>
