package com.HelixStream.HelixStreamAPI;

import org.springframework.web.bind.annotation.GetMapping;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RestController;
import org.springframework.http.ResponseEntity;

import java.io.BufferedReader;
import java.io.InputStreamReader;

@RestController
@RequestMapping("/api")
public class S3UploadController {

    @GetMapping("/s3-upload")
    public ResponseEntity<String> uploadToS3() {
        try {
            // Define the command to execute the bash script
            String[] command = {"/bin/bash", "-c", "./s3-upload.sh"};

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

            // Wait for the process to complete
            int exitCode = process.waitFor();
            if (exitCode == 0) {
                return ResponseEntity.ok("Upload successful:\n" + output.toString());
            } else {
                return ResponseEntity.status(500).body("Upload failed:\n" + output.toString());
            }
        } catch (Exception e) {
            return ResponseEntity.status(500).body("Error executing script: " + e.getMessage());
        }
    }
}