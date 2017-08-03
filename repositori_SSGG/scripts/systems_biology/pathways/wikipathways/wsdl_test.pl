use strict;
use warnings;
use SOAP::WSDL;

my %data;

my $soap = SOAP::WSDL->new(
    wsdl => 'http://www.wikipathways.org/wpi/webservice/webservice.php?wsdl',
 );

my $result = $soap->call('findPathwaysByXref', %data);

