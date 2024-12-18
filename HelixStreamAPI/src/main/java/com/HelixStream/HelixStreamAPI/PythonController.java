package com.HelixStream.HelixStreamAPI;

import org.springframework.web.bind.annotation.GetMapping;
import org.springframework.web.bind.annotation.RequestHeader;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RestController;
import org.springframework.http.ResponseEntity;

import java.io.BufferedReader;
import java.io.InputStreamReader;

@RestController
@RequestMapping("/api")
public class PythonController {

    @GetMapping("/fetch-data")
    public ResponseEntity<String> fetchData(@RequestHeader("Variant-Name") String variantName) {
        return executePythonScript("../../../../../../../../python/FetchData.py", "Fetch", variantName);
    }

    @GetMapping("/analyze-data")
    public ResponseEntity<String> analyzeData(@RequestHeader("Variant-Name") String variantName) {
        return executePythonScript("../../../../../../../../python/AnalyzeData.py", "Analyze", variantName);
    }

    @GetMapping("/parse-data")
    public ResponseEntity<String> parseData(@RequestHeader("Variant-Name") String variantName) {
        return executePythonScript("../../../../../../../../python/ParseData.py", "Parse", variantName);
    }

    private ResponseEntity<String> executePythonScript(String scriptPath, String scriptName, String variantName) {
        try {
            // Define the command to execute the Python script with the variant name as an argument
            String[] command = {"python3", scriptPath, variantName};

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
                return ResponseEntity.ok(scriptName + " successful:\n" + output.toString());
            } else {
                return ResponseEntity.status(500).body(scriptName + " failed:\n" + errorOutput.toString());
            }
        } catch (Exception e) {
            return ResponseEntity.status(500).body("Error executing " + scriptName + " script: " + e.getMessage());
        }
    }
}