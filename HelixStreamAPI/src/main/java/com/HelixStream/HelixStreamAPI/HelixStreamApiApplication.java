package com.HelixStream.HelixStreamAPI;

import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;
import org.springframework.boot.autoconfigure.jdbc.DataSourceAutoConfiguration;


@SpringBootApplication(exclude = {DataSourceAutoConfiguration.class})
public class HelixStreamApiApplication {

	public static void main(String[] args) {
		SpringApplication.run(HelixStreamApiApplication.class, args);
		System.out.println("hello");
	}

}
