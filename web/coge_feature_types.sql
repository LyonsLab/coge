SET SQL_MODE="NO_AUTO_VALUE_ON_ZERO";


/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;

--
-- Database: 'coge'
--

--
-- Dumping data for table 'feature_type'
--

INSERT INTO feature_type (feature_type_id, `name`, description) VALUES
(1, 'gene', NULL),
(2, 'mRNA', NULL),
(3, 'CDS', NULL),
(4, 'chromosome', 'entire chromosome or equivalent for genome sequencing projects'),
(5, 'tRNA', NULL),
(6, 'snoRNA', NULL),
(7, 'snRNA', NULL),
(8, 'rRNA', NULL),
(9, 'miRNA', NULL),
(10, 'miscRNA', NULL);
