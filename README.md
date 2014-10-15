The OME Framework
==

The OME Framework is a systems biology and bioinformatics framework that helps to enable a broad array of iterative and evolving systems biology workflows.

##[Installation](INSTALL.md)

###Its core features are 

1. map highthroughput experimental data to cellular components
2. map cellular components to computational models
3. map computational models to highthroughput experimental data via cellular components

![alt tag](https://raw.githubusercontent.com/steve-federowicz/om/master/docs/om_overview.png)


==

The framework is written entirely in python with great help from [SQLalchemy](http://www.sqlalchemy.org/). You can browse the schemas below to get a sense of the underlying database design or jump straight to the examples page. Enjoy!

<a href="https://www.draw.io/?url=https://raw.githubusercontent.com/steve-federowicz/om/master/docs/data_schema.xml#" target="_blank">data schema</a>

<a href="https://www.draw.io/?url=https://raw.githubusercontent.com/steve-federowicz/om/master/docs/model_schema.xml#" target="_blank">model schema</a>

<a href="https://www.draw.io/?url=https://raw.githubusercontent.com/steve-federowicz/om/master/docs/component_schema.xml#" target="_blank">component schema</a>


###Why the OME Framework?

Most biological analyses and analysis tools are focused on working with a relatively limited subset of data concerned with a specific study or question of interest.  However, as the sequencing and synthesis of DNA surges ahead it becomes clear that we will soon generate 10 replicates instead of 3, cross reference 10 historic datasets instead of 1, and most importantly begin to get a level of quantitation that will enable true biological engineering.  

The OME Framework is a small step towards a future that amasses, manages, and processes biological data at a systems scale.
