# Construct-ancestral-genome

Construct-ancestral-genome is a protocol for processing present-day genomic data, inferring local ancestry, and assembling ancestral haplotypes.

https://github.com/Shuhua-Group/Construct-ancestral-genome <br>
https://github.com/zhangxx123456/Construct-ancestral-genome <br><br><br>


### Dependencies
&#8226;&nbsp;Python (v3.6), see http://www.python.org <br>
&#8226;&nbsp;Plink (v1.9), see https://www.cog-genomics.org/plink/ <br>
&#8226;&nbsp;Plink (v2.0), see https://www.cog-genomics.org/plink/ <br>
&#8226;&nbsp;Perl (v5), see http://www.perl.org <br>
&#8226;&nbsp;ChromoPainter (v2), see https://github.com/sahwa/ChromoPainterV2 <br>
&#8226;&nbsp;BCFtools (v1.14), see https://www.htslib.org/download/ <br><br>
We create a Docker image (https://registry.hub.docker.com/r/wangbaonan/shapeit4_perl_image) that provides an environment with all software installed.<br><br><br>


### Usage
`scripts/` directory contains the core scripts used in this protocal. <br>
`Examples/` directory contains the test data and the execution scripts.  <br>
`Examples.done/` directory contains all the intermediate files and results generated after running the above pipeline, for user reference. <br>

&#8226;&nbsp;For detailed execution, please refer to the `000.run.pipeline.sh` file. <br>
&#8226;&nbsp;Before you run it, users need to modify the software paths in the `000.run.pipeline.sh` file (e.g., `python=/your/path`). <br>
&#8226;&nbsp;Then the pipeline can be quickly executed with the following command.
```
sh 000.run.pipeline.sh
```

&#8226;&nbsp;To help users run the pipeline more efficiently, we remove the data preprocessing section and the subsequent analysis for constructed ancestral samples. If users require these steps, please refer to our paper.<br><br><br>



### Citation
If you use the Construct-ancestral-genome in your research work, please cite the following papers: <br>
&nbsp;&nbsp;&nbsp;&nbsp;&#8226;&nbsp;[Protocol for reconstructing ancestral genomes from present-day samples by applying local ancestry inference](https://doi.org/10.1016/j.xpro.2024.103580) (*STAR Protocols*, 2025)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&#8226;&nbsp;[Reconstructing the ancestral gene pool to uncover the origins and genetic links of Hmong-Mien speakers](https://doi.org/10.1186/s12915-024-01838-9) (*BMC Biology*, 2024)<br><br><br>

---
By: Xiaoxi Zhang, 2024<br>
Contact: zhangxiaoxi@picb.ac.cn<br>
