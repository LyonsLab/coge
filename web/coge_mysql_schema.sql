-- MySQL dump 10.13  Distrib 5.5.49, for debian-linux-gnu (x86_64)
--
-- Host: localhost    Database: coge
-- ------------------------------------------------------
-- Server version	5.5.49-0ubuntu0.14.04.1-log

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `annotation_type`
--

DROP TABLE IF EXISTS `annotation_type`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `annotation_type` (
  `annotation_type_id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(256) NOT NULL,
  `description` varchar(1024) DEFAULT NULL,
  `annotation_type_group_id` int(11) DEFAULT NULL,
  PRIMARY KEY (`annotation_type_id`),
  KEY `name` (`name`),
  KEY `annotation_type_group_id` (`annotation_type_group_id`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `annotation_type_group`
--

DROP TABLE IF EXISTS `annotation_type_group`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `annotation_type_group` (
  `annotation_type_group_id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(256) NOT NULL,
  `description` varchar(1024) DEFAULT NULL,
  PRIMARY KEY (`annotation_type_group_id`),
  KEY `name` (`name`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `data_source`
--

DROP TABLE IF EXISTS `data_source`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `data_source` (
  `data_source_id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(256) NOT NULL,
  `description` varchar(1024) DEFAULT NULL,
  `link` text,
  PRIMARY KEY (`data_source_id`),
  KEY `name` (`name`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `dataset`
--

DROP TABLE IF EXISTS `dataset`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `dataset` (
  `dataset_id` int(11) NOT NULL AUTO_INCREMENT,
  `data_source_id` int(11) NOT NULL DEFAULT '0',
  `name` varchar(255) NOT NULL,
  `description` varchar(1024) DEFAULT NULL,
  `version` varchar(50) DEFAULT NULL,
  `link` text,
  `creator_id` int(11) NOT NULL,
  `date` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  `restricted` int(1) NOT NULL DEFAULT '0',
  `deleted` tinyint(1) NOT NULL DEFAULT '0',
  PRIMARY KEY (`dataset_id`),
  KEY `data_source_id` (`data_source_id`),
  KEY `name` (`name`),
  KEY `restricted` (`restricted`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `dataset_connector`
--

DROP TABLE IF EXISTS `dataset_connector`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `dataset_connector` (
  `dataset_connector_id` int(11) NOT NULL AUTO_INCREMENT,
  `dataset_id` int(11) NOT NULL,
  `genome_id` int(11) NOT NULL,
  PRIMARY KEY (`dataset_connector_id`),
  KEY `dataset_id` (`dataset_id`),
  KEY `dataset_id_2` (`dataset_id`,`genome_id`),
  KEY `dataset_group_id` (`genome_id`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `experiment`
--

DROP TABLE IF EXISTS `experiment`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `experiment` (
  `experiment_id` int(11) NOT NULL AUTO_INCREMENT,
  `genome_id` int(11) NOT NULL,
  `data_source_id` int(11) NOT NULL,
  `data_type` tinyint(1) DEFAULT NULL,
  `name` varchar(255) NOT NULL,
  `description` text,
  `version` varchar(50) NOT NULL,
  `storage_path` varchar(255) NOT NULL,
  `restricted` int(1) NOT NULL DEFAULT '0',
  `row_count` int(10) DEFAULT NULL,
  `creator_id` int(11) NOT NULL,
  `date` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  `deleted` tinyint(1) NOT NULL DEFAULT '0',
  PRIMARY KEY (`experiment_id`),
  KEY `dataset_group_id` (`genome_id`),
  KEY `data_source_id` (`data_source_id`),
  KEY `restricted` (`restricted`),
  FULLTEXT KEY `name` (`name`),
  FULLTEXT KEY `description` (`description`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `experiment_annotation`
--

DROP TABLE IF EXISTS `experiment_annotation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `experiment_annotation` (
  `experiment_annotation_id` int(11) NOT NULL AUTO_INCREMENT,
  `experiment_id` int(11) NOT NULL,
  `annotation_type_id` int(11) NOT NULL,
  `annotation` text NOT NULL,
  `link` text,
  `image_id` int(11) DEFAULT NULL,
  `locked` int(1) NOT NULL DEFAULT '0',
  PRIMARY KEY (`experiment_annotation_id`),
  KEY `experiment_id` (`experiment_id`),
  KEY `annotation_type_id` (`annotation_type_id`),
  FULLTEXT KEY `annotation` (`annotation`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `experiment_type`
--

DROP TABLE IF EXISTS `experiment_type`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `experiment_type` (
  `experiment_type_id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(255) NOT NULL,
  PRIMARY KEY (`experiment_type_id`),
  FULLTEXT KEY `name` (`name`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `experiment_type_connector`
--

DROP TABLE IF EXISTS `experiment_type_connector`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `experiment_type_connector` (
  `experiment_type_connector_id` int(11) NOT NULL AUTO_INCREMENT,
  `experiment_type_id` int(11) NOT NULL,
  `experiment_id` int(11) NOT NULL,
  PRIMARY KEY (`experiment_type_connector_id`),
  KEY `experiment_type_id` (`experiment_type_id`),
  KEY `experiment_id` (`experiment_id`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `feature`
--

DROP TABLE IF EXISTS `feature`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `feature` (
  `feature_id` int(11) NOT NULL AUTO_INCREMENT,
  `feature_type_id` int(11) NOT NULL,
  `dataset_id` int(11) NOT NULL,
  `start` int(11) DEFAULT NULL,
  `stop` int(11) DEFAULT NULL,
  `strand` tinyint(4) DEFAULT NULL,
  `chromosome` varchar(255) DEFAULT NULL,
  `access_count` int(10) NOT NULL DEFAULT '0',
  PRIMARY KEY (`feature_id`),
  KEY `feature_type_id` (`feature_type_id`),
  KEY `dataset_id` (`dataset_id`),
  KEY `start` (`start`),
  KEY `stop` (`stop`),
  KEY `chromosome` (`chromosome`),
  KEY `dataset_id_2` (`dataset_id`,`chromosome`),
  KEY `dataset_feature_type` (`feature_type_id`,`dataset_id`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8
/*!50100 PARTITION BY HASH (feature_id)
PARTITIONS 101 */;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `feature_annotation`
--

DROP TABLE IF EXISTS `feature_annotation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `feature_annotation` (
  `feature_annotation_id` int(11) NOT NULL AUTO_INCREMENT,
  `annotation` text,
  `feature_id` int(11) NOT NULL DEFAULT '0',
  `annotation_type_id` int(11) NOT NULL DEFAULT '0',
  `link` varchar(1024) DEFAULT NULL,
  PRIMARY KEY (`feature_annotation_id`),
  KEY `feature_id` (`feature_id`),
  KEY `annotation_type_id` (`annotation_type_id`),
  FULLTEXT KEY `annotation` (`annotation`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `feature_name`
--

DROP TABLE IF EXISTS `feature_name`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `feature_name` (
  `feature_name_id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(255) NOT NULL,
  `description` varchar(1024) DEFAULT NULL,
  `feature_id` int(10) NOT NULL DEFAULT '0',
  `primary_name` tinyint(1) NOT NULL DEFAULT '0',
  PRIMARY KEY (`feature_name_id`),
  KEY `name` (`name`),
  KEY `feature_id` (`feature_id`),
  FULLTEXT KEY `name_2` (`name`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `feature_type`
--

DROP TABLE IF EXISTS `feature_type`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `feature_type` (
  `feature_type_id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(255) NOT NULL,
  `description` varchar(1024) DEFAULT NULL,
  PRIMARY KEY (`feature_type_id`),
  KEY `name` (`name`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `genome`
--

DROP TABLE IF EXISTS `genome`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `genome` (
  `genome_id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(200) DEFAULT NULL,
  `description` varchar(255) DEFAULT NULL,
  `version` varchar(50) NOT NULL,
  `organism_id` int(11) NOT NULL,
  `genomic_sequence_type_id` int(11) NOT NULL,
  `file_path` varchar(255) DEFAULT NULL,
  `restricted` tinyint(1) NOT NULL DEFAULT '0',
  `access_count` int(10) DEFAULT NULL,
  `message` text,
  `link` text,
  `deleted` tinyint(1) NOT NULL DEFAULT '0',
  `creator_id` int(11) NOT NULL,
  `date` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (`genome_id`),
  KEY `name` (`name`),
  KEY `organism_id` (`organism_id`),
  KEY `genome_sequence_type_id` (`genomic_sequence_type_id`),
  KEY `version` (`version`,`organism_id`),
  KEY `restricted` (`restricted`),
  KEY `creator_id` (`creator_id`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `genome_annotation`
--

DROP TABLE IF EXISTS `genome_annotation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `genome_annotation` (
  `genome_annotation_id` int(11) NOT NULL AUTO_INCREMENT,
  `genome_id` int(11) NOT NULL,
  `annotation_type_id` int(11) NOT NULL,
  `annotation` text NOT NULL,
  `link` text,
  `image_id` int(11) DEFAULT NULL,
  `locked` int(1) NOT NULL DEFAULT '0',
  PRIMARY KEY (`genome_annotation_id`),
  KEY `genome_id` (`genome_id`),
  KEY `annotation_type_id` (`annotation_type_id`),
  FULLTEXT KEY `annotation` (`annotation`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `genomic_sequence`
--

DROP TABLE IF EXISTS `genomic_sequence`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `genomic_sequence` (
  `genomic_sequence_id` int(11) NOT NULL AUTO_INCREMENT,
  `sequence_length` int(11) NOT NULL,
  `chromosome` varchar(255) NOT NULL,
  `genome_id` int(11) NOT NULL,
  PRIMARY KEY (`genomic_sequence_id`),
  KEY `chromosome` (`chromosome`),
  KEY `dataset_group_id` (`genome_id`),
  KEY `chromosome_2` (`chromosome`,`genome_id`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `genomic_sequence_type`
--

DROP TABLE IF EXISTS `genomic_sequence_type`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `genomic_sequence_type` (
  `genomic_sequence_type_id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(256) NOT NULL,
  `description` varchar(1024) DEFAULT NULL,
  PRIMARY KEY (`genomic_sequence_type_id`),
  KEY `name` (`name`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `image`
--

DROP TABLE IF EXISTS `image`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `image` (
  `image_id` int(11) NOT NULL AUTO_INCREMENT,
  `filename` varchar(255) NOT NULL,
  `image` longblob NOT NULL,
  PRIMARY KEY (`image_id`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `job`
--

DROP TABLE IF EXISTS `job`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `job` (
  `job_id` int(11) NOT NULL AUTO_INCREMENT,
  `user_id` int(11) NOT NULL,
  `page` varchar(255) NOT NULL,
  `start_time` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  `end_time` timestamp NULL DEFAULT NULL,
  `link` varchar(255) DEFAULT NULL,
  `status` tinyint(5) NOT NULL DEFAULT '0',
  `type` tinyint(5) NOT NULL DEFAULT '0',
  `process_id` int(5) NOT NULL,
  `log_id` int(11) DEFAULT NULL,
  PRIMARY KEY (`job_id`),
  KEY `user_id` (`user_id`),
  FULLTEXT KEY `page` (`page`),
  FULLTEXT KEY `link` (`link`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `list`
--

DROP TABLE IF EXISTS `list`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `list` (
  `list_id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(255) NOT NULL,
  `description` varchar(1024) DEFAULT NULL,
  `list_type_id` int(11) NOT NULL,
  `creator_id` int(11) NOT NULL,
  `date` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  `restricted` tinyint(1) NOT NULL DEFAULT '0',
  `locked` int(1) NOT NULL DEFAULT '0',
  `deleted` tinyint(1) NOT NULL DEFAULT '0',
  PRIMARY KEY (`list_id`),
  KEY `name` (`name`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `list_annotation`
--

DROP TABLE IF EXISTS `list_annotation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `list_annotation` (
  `list_annotation_id` int(11) NOT NULL AUTO_INCREMENT,
  `list_id` int(11) NOT NULL,
  `annotation_type_id` int(11) NOT NULL,
  `annotation` text NOT NULL,
  `link` text,
  `image_id` int(11) DEFAULT NULL,
  `locked` int(1) NOT NULL DEFAULT '0',
  PRIMARY KEY (`list_annotation_id`),
  KEY `list_id` (`list_id`),
  KEY `annotation_type_id` (`annotation_type_id`),
  FULLTEXT KEY `annotation` (`annotation`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `list_connector`
--

DROP TABLE IF EXISTS `list_connector`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `list_connector` (
  `list_connector_id` int(11) NOT NULL AUTO_INCREMENT,
  `parent_id` int(11) DEFAULT NULL,
  `child_id` int(11) NOT NULL,
  `child_type` tinyint(1) NOT NULL,
  PRIMARY KEY (`list_connector_id`),
  KEY `parent_id` (`parent_id`,`child_id`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `list_type`
--

DROP TABLE IF EXISTS `list_type`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `list_type` (
  `list_type_id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(255) NOT NULL,
  `description` varchar(1024) DEFAULT NULL,
  PRIMARY KEY (`list_type_id`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `location`
--

DROP TABLE IF EXISTS `location`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `location` (
  `location_id` int(10) NOT NULL AUTO_INCREMENT,
  `start` int(10) NOT NULL DEFAULT '0',
  `stop` int(10) NOT NULL DEFAULT '0',
  `chromosome` varchar(255) NOT NULL,
  `feature_id` int(10) NOT NULL DEFAULT '0',
  `strand` tinyint(4) NOT NULL,
  PRIMARY KEY (`location_id`),
  KEY `feature_id` (`feature_id`),
  KEY `chromosome` (`chromosome`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8
/*!50100 PARTITION BY HASH (location_id)
PARTITIONS 101 */;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `log`
--

DROP TABLE IF EXISTS `log`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `log` (
  `log_id` int(11) NOT NULL AUTO_INCREMENT,
  `time` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  `user_id` int(11) NOT NULL DEFAULT '0',
  `type` tinyint(1) NOT NULL DEFAULT '0',
  `page` varchar(255) NOT NULL,
  `description` varchar(255) DEFAULT NULL,
  `link` varchar(255) DEFAULT NULL,
  `status` int(1) DEFAULT '0',
  `comment` varchar(255) DEFAULT NULL,
  `parent_id` int(11) DEFAULT NULL,
  `parent_type` tinyint(4) DEFAULT NULL,
  PRIMARY KEY (`log_id`),
  KEY `user_id` (`user_id`),
  KEY `type` (`type`),
  KEY `user_id_2` (`user_id`,`type`),
  KEY `time` (`time`,`user_id`,`type`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `man_proof_cns_coge`
--

DROP TABLE IF EXISTS `man_proof_cns_coge`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `man_proof_cns_coge` (
  `feature_id` int(10) NOT NULL DEFAULT '0',
  `feature_type_id` int(10) NOT NULL DEFAULT '0',
  `dataset_id` int(10) NOT NULL DEFAULT '0',
  `start` int(11) DEFAULT NULL,
  `stop` int(11) DEFAULT NULL,
  `strand` tinyint(4) DEFAULT NULL,
  `chromosome` varchar(255) DEFAULT NULL,
  `access_count` int(10) NOT NULL DEFAULT '0',
  `name` varchar(255) NOT NULL,
  `annotation` text
) ENGINE=MyISAM DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `organism`
--

DROP TABLE IF EXISTS `organism`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `organism` (
  `organism_id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(255) NOT NULL,
  `description` varchar(1024) DEFAULT NULL,
  `restricted` int(1) NOT NULL DEFAULT '0',
  `access_count` int(10) NOT NULL DEFAULT '0',
  PRIMARY KEY (`organism_id`),
  KEY `name` (`name`),
  KEY `restricted` (`restricted`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `permission`
--

DROP TABLE IF EXISTS `permission`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `permission` (
  `permission_id` int(10) NOT NULL AUTO_INCREMENT,
  `name` varchar(256) NOT NULL,
  `description` varchar(1024) DEFAULT NULL,
  PRIMARY KEY (`permission_id`),
  KEY `name` (`name`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `role`
--

DROP TABLE IF EXISTS `role`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `role` (
  `role_id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(255) NOT NULL,
  `description` longtext,
  PRIMARY KEY (`role_id`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `role_permission_connector`
--

DROP TABLE IF EXISTS `role_permission_connector`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `role_permission_connector` (
  `role_permission_connector_id` int(11) NOT NULL AUTO_INCREMENT,
  `role_id` int(11) NOT NULL,
  `permission_id` int(11) NOT NULL,
  PRIMARY KEY (`role_permission_connector_id`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `user`
--

DROP TABLE IF EXISTS `user`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `user` (
  `user_id` int(11) NOT NULL AUTO_INCREMENT,
  `user_name` varchar(50) NOT NULL,
  `first_name` varchar(50) NOT NULL,
  `last_name` varchar(50) NOT NULL,
  `email` varchar(50) DEFAULT NULL,
  `description` varchar(255) DEFAULT NULL,
  `image_id` int(11) DEFAULT NULL,
  `date` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  `admin` tinyint(1) NOT NULL DEFAULT '0',
  PRIMARY KEY (`user_id`),
  KEY `user_name` (`user_name`),
  KEY `first_name` (`first_name`),
  KEY `last_name` (`last_name`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `user_connector`
--

DROP TABLE IF EXISTS `user_connector`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `user_connector` (
  `user_connector_id` int(11) NOT NULL AUTO_INCREMENT,
  `parent_id` int(11) DEFAULT NULL,
  `parent_type` tinyint(1) NOT NULL,
  `child_id` int(11) NOT NULL,
  `child_type` tinyint(1) NOT NULL,
  `role_id` int(11) NOT NULL DEFAULT '4',
  PRIMARY KEY (`user_connector_id`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `user_group`
--

DROP TABLE IF EXISTS `user_group`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `user_group` (
  `user_group_id` int(10) NOT NULL AUTO_INCREMENT,
  `creator_user_id` int(11) NOT NULL,
  `date` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  `name` varchar(255) NOT NULL,
  `description` varchar(1024) NOT NULL,
  `role_id` int(11) NOT NULL DEFAULT '4',
  `locked` int(1) NOT NULL DEFAULT '0',
  `deleted` tinyint(1) NOT NULL DEFAULT '0',
  PRIMARY KEY (`user_group_id`),
  KEY `name` (`name`),
  KEY `role_id` (`role_id`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `user_session`
--

DROP TABLE IF EXISTS `user_session`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `user_session` (
  `user_session_id` int(11) NOT NULL AUTO_INCREMENT,
  `user_id` int(11) NOT NULL,
  `date` datetime NOT NULL,
  `session` varchar(22) NOT NULL,
  PRIMARY KEY (`user_session_id`),
  KEY `user_id` (`user_id`),
  KEY `session` (`session`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `web_preferences`
--

DROP TABLE IF EXISTS `web_preferences`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `web_preferences` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `user_id` int(11) NOT NULL,
  `page` varchar(255) NOT NULL,
  `options` text NOT NULL,
  PRIMARY KEY (`id`),
  KEY `user_id` (`user_id`),
  KEY `page` (`page`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `work`
--

DROP TABLE IF EXISTS `work`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `work` (
  `work_id` int(11) NOT NULL AUTO_INCREMENT,
  `user_id` int(11) NOT NULL,
  `page` varchar(255) NOT NULL,
  `parameter` text,
  `name` varchar(255) DEFAULT NULL,
  `description` varchar(1024) DEFAULT NULL,
  `note` text,
  `archive` smallint(1) DEFAULT '0' COMMENT '1 for archive',
  `image_id` int(11) DEFAULT NULL,
  `date` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  `link` text,
  PRIMARY KEY (`work_id`),
  KEY `user_id` (`user_id`),
  KEY `page` (`page`),
  KEY `image_id` (`image_id`),
  FULLTEXT KEY `notes` (`note`),
  FULLTEXT KEY `description` (`description`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `work_order`
--

DROP TABLE IF EXISTS `work_order`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `work_order` (
  `work_order_id` int(11) NOT NULL AUTO_INCREMENT,
  `workflow_id` int(11) NOT NULL,
  `work_id` int(11) NOT NULL,
  `work_order` int(4) NOT NULL DEFAULT '1' COMMENT 'order in the workflow',
  `date` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`work_order_id`),
  KEY `workflow_id` (`workflow_id`),
  KEY `work_id` (`work_id`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `workflow`
--

DROP TABLE IF EXISTS `workflow`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `workflow` (
  `workflow_id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(255) NOT NULL,
  `description` varchar(1024) DEFAULT NULL,
  `user_id` int(11) NOT NULL,
  `date` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  `note` text,
  `link` text,
  `image_id` int(11) DEFAULT NULL,
  PRIMARY KEY (`workflow_id`),
  KEY `user_id` (`user_id`),
  KEY `name` (`name`),
  FULLTEXT KEY `description` (`description`)
) ENGINE=MyISAM  DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2016-06-03 10:16:36
