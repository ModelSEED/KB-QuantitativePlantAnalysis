/*
A KBase module: QuantitativePlantAnalysis
*/

module QuantitativePlantAnalysis {
    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;

    /*
        Compute plant biomass yield from glucose
    */
    funcdef compute_plant_biomass_yield(mapping<string,UnspecifiedObject> params) returns (ReportResults output) authentication required;

};
