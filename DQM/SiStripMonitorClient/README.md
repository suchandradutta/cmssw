###Macros and Scripts
####moduleOccupancyPlots
`moduleOccupancyPlots.sh <Datataking_period> <Dataset_type> <runnumber> <modulelistfile> <user certificate file> <user key file>`

This script produces a root file and png files of the occupancy plots of the modules selected with the file 
`<modulelistfile>` which contains a list of detid from run `<runnumber>` as found in the official DQM file. The 
`<Datataking_period>` and the `<Dataset_type>` have to match the names used in this site: 
https://cmsweb.cern.ch/dqm/offline/data/browse/ROOT/OfflineData/

To access the DQM file a valid certificate and key has to be provided to the script

- New functionality added to create a summary histogram adding individual module histograms. The corresponding png file
  is also created and the summary histogram is added in the output root file. The Summary histogram is scaled by the # of 
  events.

- Also added option (=8) to create re-ordered (as of APV sequence) occupancy histograms from module occupancy histograms
  using following formula where ival is the sequence number and chan is corresponding APV channel number

    int chan=int(32*fmod(int(fmod(ival,256.)/2.),4.) +
              8*int(int(fmod(ival,256.)/2.)/4.) -
              31*int(int(fmod(ival,256.)/2.)/16.) +
              fmod(fmod(ival,256.),2.)*128 +
              int(ival/256)*256);
