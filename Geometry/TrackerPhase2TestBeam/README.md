                      Tracker Geometry                                      ||                 PixelTelescope geometry                                    ||       Notes
                                                                            ||                                                                            ||
                                                                            ||                                                                            ||
                                                                            ||                                                                            ||
                                                                            ||                                                                            ||
Geometry.CMSCommonData.cmsExtendedGeometry2023D21XML_cfi                    ||     Geometry.TrackerPhase2TestBeam.Phase2TestBeamGeometryXML_cfi           ||  DD Geometry.
                                                                            ||                                                                            ||                          
Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cfi               ||     Geometry.TrackerPhase2TestBeam.telescopeGeometryNumbering_cfi          ||  Geometry is ordered and DetIds are assigned.
                                                                            ||                                                                            ||
Geometry.TrackerGeometryBuilder.trackerParameters_cfi                       ||     Geometry.TrackerPhase2TestBeam.telescopeParameters_cfi                 ||  Read parameters from DD (actually, DetId sheme only).
                                                                            ||                                                                            ||
Geometry.TrackerNumberingBuilder.trackerTopology_cfi                        ||     Geometry.TrackerPhase2TestBeam.telescopeTopology_cfi                   ||  Allow to get the layer, or plane, or whether a sensor is inner or outer, etc.. from a given DetId.
                                                                            ||                                                                            ||
Geometry.TrackerGeometryBuilder.trackerGeometry_cfi                         ||     Geometry.TrackerGeometryBuilder.trackerGeometry_cfi                    ||  Full geometry, as used by the Digitizer. 
                                                                            ||                                                                            ||  NB: Already operational on telescope geometry, but needs to be renamed to telescopeGeometry_cfi and moved to Phase2TestBeam package.
                                                                            ||                                                                            ||
SLHCUpgradeSimulations.Geometry.fakeConditions_phase2TkTiltedBase_cff       ||                          ?                                                 ||  Not tuned yet.



The work so far more or less create the equivalent of the full TrackerNumberingBuilder and TrackerGeometry packages for the telescope.
