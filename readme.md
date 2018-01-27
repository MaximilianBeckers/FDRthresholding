# False discovery rate control of cryo-EM maps

Confidence maps are complementary maps generated from cryo-EM maps by means of statistical hypothesis testing and subsequent FDR control. They allow thresholding of EM maps based on the expected amount of background noise visible at the respective threshold and thus allow rigorous error assesment of visible features in the density. 
Additional postprocessing like local filtering or local amplitude scaling (LocScale) can be incorporated in the framework in order to increase the power.


## Getting Started


### Prerequisites

The software is written in Python and itself dependent on EMAN2, NumPy and SciPy libraries. 

For incorporation of the local ampltidue scaling (LocScale) procedue, additional LocScale libraries are needed, as described here:  https://git.embl.de/jakobi/LocScale/wikis/home


### Installing

The software consists of four scripts, FDRcontrol.py, mapUtil.py, FDRUtil.py and locscaleUtil.py, that just have to be copied to your computer.

Once you have a working EMAN2 installation ( http://blake.bcm.edu/emanwiki/EMAN2 ), the software can be simply run by using your Python version that comes with EMAN2:   

For example if /programs/x86_64-linux/eman2/2.2/bin/python is the respective EMAN2 python installation: 
```
/programs/x86_64-linux/eman2/2.2/bin/python FDRcontrol.py -em yourMap.mrc -p thePixelSize
```

If you want to use the LocScale feature, you should have a running LocScale version, together with MPI,  as described here:  https://git.embl.de/jakobi/LocScale/wikis/home.  


## How to use






### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone who's code was used
* Inspiration
* etc
