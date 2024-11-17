package com.HelixStream.HelixStreamAPI;

import org.springframework.web.bind.annotation.GetMapping;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RestController;
import org.springframework.http.ResponseEntity;

import java.io.BufferedReader;
import java.io.InputStreamReader;

@RestController
@RequestMapping("/api")
public class FetchController {

    @GetMapping("/fetch-data")
    public ResponseEntity<String> fetchData() {
        try {
            // Define the command to execute the Python script
            String[] command = {"python3", "../../../../../../../../python/FetchData.py"};

            // Create a ProcessBuilder to execute the command
            ProcessBuilder processBuilder = new ProcessBuilder(command);

            // Start the process
            Process process = processBuilder.start();

            // Capture the output of the script
            BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
            StringBuilder output = new StringBuilder();
            String line;
            while ((line = reader.readLine()) != null) {
                output.append(line).append("\n");
            }

            // Capture any errors from the script
            BufferedReader errorReader = new BufferedReader(new InputStreamReader(process.getErrorStream()));
            StringBuilder errorOutput = new StringBuilder();
            while ((line = errorReader.readLine()) != null) {
                errorOutput.append(line).append("\n");
            }

            // Wait for the process to complete
            int exitCode = process.waitFor();
            if (exitCode == 0) {
                return ResponseEntity.ok("Fetch successful:\n" + output.toString());
            } else {
                return ResponseEntity.status(500).body("Fetch failed:\n" + errorOutput.toString());
            }
        } catch (Exception e) {
            return ResponseEntity.status(500).body("Error executing script: " + e.getMessage());
        }
    }
}