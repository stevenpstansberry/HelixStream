package com.HelixStream.HelixStreamAPI;

import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;
import org.springframework.boot.autoconfigure.jdbc.DataSourceAutoConfiguration;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

@SpringBootApplication(exclude = {DataSourceAutoConfiguration.class})
public class HelixStreamApiApplication {

	private static final Logger logger = LoggerFactory.getLogger(HelixStreamApiApplication.class);

	public static void main(String[] args) {
		SpringApplication.run(HelixStreamApiApplication.class, args);

		logger.info("=======================================================");
		logger.info("     HelixStream API Server Started Successfully!      ");
		logger.info("=======================================================");
		logger.info("Application Name: HelixStream API");
		logger.info("Version: 1.0.0");
		logger.info("Environment: {}", System.getProperty("spring.profiles.active", "default"));
		logger.info("Listening on port: {}", System.getProperty("server.port", "8080"));
		logger.info("Start Time: {}", java.time.LocalDateTime.now());
		logger.info("=======================================================");
	}
}
