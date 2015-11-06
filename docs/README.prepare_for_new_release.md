----------------------------------------------------------------------------
Preparing ALLENMINER release v2.1

Oct 28, 2015.
Fred P. Davis
----------------------------------------------------------------------------


o Preparing release 2.1
-----------------------

* code change: new atlas URLs

* code change: new Atlas ontology file formats

* code change: new developmental atlas gene list web page URL and formats

    * OLD url: http://developingmouse.brain-map.org/data/search/gene/index.xml?term=*

    * NEW url: http://developingmouse.brain-map.org/api/v2/data/Gene/query.xml?query=*0*&startRow=1300

    * NOT TRUE: http://developingmouse.brain-map.org/api/v2/data/Gene/query.xml will give all genes

    * have to query with individual terms like doing for adult atlas

    * HOLY CRAP, complicated API now, both for adult and devel atlases:
    * http://help.brain-map.org/display/mousebrain/API
    * http://help.brain-map.org/display/devmouse/API
    * http://help.brain-map.org/display/api/Connected+Services+and+Pipes#ConnectedServicesandPipes-service%3A%3Amousecorrelation

* download atlases and expression data

    * adult: done
    * spinal: done
    * devel: 

* build fastsearch index
   perl allenminer.pl -mode fastsearch -prep_data 1

* post fastsearch files on Zenodo

* load to github
