#include "Record.h"

#include <iostream>

#include "utilities.h"
#include "classify.h"
#include "json.hpp"
#include "artefact.h"
#include "strand.h"
#include <iomanip> // For formatting the output table

/*
    initialises an empty record
*/
Record::Record()
{

}

/*
    initialises a record
*/
Record::Record
(
    std::string identifier,
    std::string sequence,
    std::string quality
)
{
    // std::cerr << "created record:\n"
    //         << "\tidentifier: " << identifier << "\n"
    //         << "\tsequence: " << sequence << "\n" 
    //         << "\tquality: " << quality << "\n"
            // << "\tmethod: " << method << "\n";

    this->identifier = identifier;
    this->sequence = sequence;
    this->quality = quality;
    this->artefact = artefact::none;
}

/*
    classifies the record using a given pipeline
*/
void
Record::classify(Pipeline& pipeline)
{
    auto result = classifyPipeline(this->sequence, pipeline);
    this->strand = result.strand;
    this->artefact = result.artefact;
}

/*
    prints out a record with correct formatting
*/
std::string
Record::printFq()
{
    // first, write down if it's an artefact
    if (this->artefact != artefact::none) {
        this->identifier += " artefact=" + artefact::getName(this->artefact) + " ";
    }

    // first, add on the strand tag
    this->identifier += " strand=";
    this->identifier.push_back(this->strand);
    
    // if it's a reverse read, we need to take the reverse of some things
    if (this->strand == '-') {
        this->sequence = reverseComplement(this->sequence);
        this->quality = reverse(this->quality);
    }
    
    // then construct the formatted record
    return this->identifier + "\n" + this->sequence + "\n+\n" + this->quality + "\n";
}

/*
    This method reports a breakdown of strand and artefact relationships in a matrix format.
    It requires a collection of records to summarize their properties.
*/
void ReportMatrix(const std::vector<Record>& records) {
    // Define strand and artefact categories
    const std::vector<std::string> strands = {"+", "-", "?"};
    const std::vector<std::string> artefacts = {"no artefact", "RTP-RTP", "TSO-TSO"};

    // Matrix-like data structure to count occurrences
    std::unordered_map<std::string, std::unordered_map<std::string, int>> matrix;

    // Initialize matrix with zeros
    for (const auto& strand : strands) {
        for (const auto& artefact : artefacts) {
            matrix[strand][artefact] = 0;
        }
    }

    // Populate matrix with counts from records
for (const auto& record : records) {
        std::string strand(1, record.strand);
        std::string artefact = artefact::getName(record.artefact);
        matrix[strand][artefact]++;
    }

    // Print the matrix
    std::cout << std::setw(12) << "Strand\\Artefact";
    for (const auto& artefact : artefacts) {
        std::cout << std::setw(12) << artefact;
    }
    std::cout << std::endl;

    for (const auto& strand : strands) {
        std::cout << std::setw(12) << strand;
        for (const auto& artefact : artefacts) {
            std::cout << std::setw(12) << matrix[strand][artefact];
        }
        std::cout << std::endl;
    }

    // Additional totals or relationships
    int totalRecords = 0;
    for (const auto& strandEntry : matrix) {
        for (const auto& count : strandEntry.second) {
            totalRecords += count.second;
        }
    }
    std::cout << "Total Records: " << totalRecords << std::endl;
}
